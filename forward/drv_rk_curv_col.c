/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"
#include "drv_rk_curv_col.h"
#include "sv_curv_col_el.h"
#include "sv_curv_col_el_iso.h"
#include "sv_curv_col_el_vti.h"
#include "sv_curv_col_el_aniso.h"
#include "sv_curv_col_ac_iso.h"

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 ******************************************************************************/

int
drv_rk_curv_col_allstep(
  fd_t        *fd,
  gd_t        *gd,
  gdcurv_metric_t *metric,
  md_t      *md,
  src_t      *src,
  bdry_t *bdry,
  wav_t  *wav,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose)
{
  // retrieve from struct
  int num_rk_stages = fd->num_rk_stages;
  float *rk_a = fd->rk_a;
  float *rk_b = fd->rk_b;

  int num_of_pairs     = fd->num_of_pairs;
  int fdx_max_half_len = fd->fdx_max_half_len;
  int fdz_max_len      = fd->fdz_max_len;
  int num_of_fdz_op    = fd->num_of_fdz_op;

  // local allocated array
  char ou_file[CONST_MAX_STRLEN];

  // local pointer
  float *restrict w_cur;
  float *restrict w_pre;
  float *restrict w_rhs;
  float *restrict w_end;
  float *restrict w_tmp;

  int   ipair, istage;
  float t_cur;
  float t_end; // time after this loop for nc output
  // for mpi message

  // create snapshot nc output files
  if (verbose>0) fprintf(stdout,"prepare snap nc output ...\n"); 
  iosnap_nc_t  iosnap_nc;
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    io_snap_nc_create_ac(iosnap, &iosnap_nc);
  } else {
    io_snap_nc_create(iosnap, &iosnap_nc);
  }

  // get wavefield
  w_pre = wav->v4d + wav->siz_ilevel * 0; // previous level at n
  w_tmp = wav->v4d + wav->siz_ilevel * 1; // intermidate value
  w_rhs = wav->v4d + wav->siz_ilevel * 2; // for rhs
  w_end = wav->v4d + wav->siz_ilevel * 3; // end level at n+1

  // set pml for rk
  if(bdry->is_enable_pml == 1)
  {
    for (int idim=0; idim<CONST_NDIM; idim++) {
      for (int iside=0; iside<2; iside++) {
        if (bdry->is_sides_pml[idim][iside]==1) {
          bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
          auxvar->pre = auxvar->var + auxvar->siz_ilevel * 0;
          auxvar->tmp = auxvar->var + auxvar->siz_ilevel * 1;
          auxvar->rhs = auxvar->var + auxvar->siz_ilevel * 2;
          auxvar->end = auxvar->var + auxvar->siz_ilevel * 3;
        }
      }
    }
  }

  // calculate conversion matrix for free surface
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    if(md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
    {
      sv_curv_col_el_iso_dvh2dvz(gd,metric,md,bdry,verbose);
    } 
    else if(md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) 
    {
      sv_curv_col_el_aniso_dvh2dvz(gd,metric,md,bdry,verbose);
    }
    if(md->medium_type == CONST_MEDIUM_ELASTIC_VTI)
    {
      sv_curv_col_el_vti_dvh2dvz(gd,metric,md,bdry,verbose);
    } 
    else if(md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) 
    {
      // no need
    }
  }

  //--------------------------------------------------------
  // time loop
  //--------------------------------------------------------

  if (verbose>0) fprintf(stdout,"start time loop ...\n"); 

  for (int it=0; it<nt_total; it++)
  {
    t_cur  = it * dt + t0;
    t_end = t_cur +dt;

    if (verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // mod to get ipair
    ipair = it % num_of_pairs;
    if (verbose>10) fprintf(stdout, " --> ipair=%d\n",ipair);

    // loop RK stages for one step
    for (istage=0; istage<num_rk_stages; istage++)
    {
      if (verbose>10) fprintf(stdout, " --> istage=%d\n",istage);

      // use pointer to avoid 1 copy for previous level value
      if (istage==0) {
        w_cur = w_pre;
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              bdry->auxvar[idim][iside].cur = bdry->auxvar[idim][iside].pre;
            }
          }
        }
      }
      else
      {
        w_cur = w_tmp;
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              bdry->auxvar[idim][iside].cur = bdry->auxvar[idim][iside].tmp;
            }
          }
        }
      }

      // set src_t time
      src_set_time(src, it, istage);

      // compute rhs
      if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
      {
          sv_curv_col_el_iso_onestage(
              w_cur,w_rhs,wav,
              gd, metric, md, bdry, src,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI)
      {
          sv_curv_col_el_vti_onestage(
              w_cur,w_rhs,wav,
              gd, metric, md, bdry, src,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO)
      {
          sv_curv_col_el_aniso_onestage(
              w_cur,w_rhs,wav,
              gd, metric, md, bdry, src,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
      } else if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO)
      {
          sv_curv_col_ac_iso_onestage(
              w_cur,w_rhs,wav,
              gd, metric, md, bdry, src,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
      }

      // rk start
      if (istage==0)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wav->siz_ilevel; iptr++) {
            w_tmp[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }

        // apply Qs
        //if (md->visco_type == CONST_VISCO_GRAVES) {
        //  sv_curv_graves_Qs(w_tmp, wave->ncmp, gd, md);
        //}

        // pml_tmp
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry->is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
                for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                  auxvar->tmp[iptr] = auxvar->pre[iptr] + coef_a * auxvar->rhs[iptr];
                }
              }
            }
          }
        }

        // w_end
        for (size_t iptr=0; iptr < wav->siz_ilevel; iptr++) {
            w_end[iptr] = w_pre[iptr] + coef_b * w_rhs[iptr];
        }
        // pml_end
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry->is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
                for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                  auxvar->end[iptr] = auxvar->pre[iptr] + coef_b * auxvar->rhs[iptr];
                }
              }
            }
          }
        }
      }
      else if (istage<num_rk_stages-1)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wav->siz_ilevel; iptr++) {
            w_tmp[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }

        // apply Qs
        //if (md->visco_type == CONST_VISCO_GRAVES) {
        //  sv_curv_graves_Qs(w_tmp, wave->ncmp, gd, md);
        //}

        // pml_tmp
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry->is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
                for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                  auxvar->tmp[iptr] = auxvar->pre[iptr] + coef_a * auxvar->rhs[iptr];
                }
              }
            }
          }
        }

        // w_end
        for (size_t iptr=0; iptr < wav->siz_ilevel; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        // pml_end
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry->is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
                for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                  auxvar->end[iptr] += coef_b * auxvar->rhs[iptr];
                }
              }
            }
          }
        }
      }
      else // last stage
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wav->siz_ilevel; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }

        // apply Qs
        if (md->visco_type == CONST_VISCO_GRAVES_QS) {
          sv_curv_graves_Qs(w_end, wav->ncmp, dt, gd, md);
        }
        
        // pml_end
        if(bdry->is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry->is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
                for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                  auxvar->end[iptr] += coef_b * auxvar->rhs[iptr];
                }
              }
            }
          }
        }
      }

    } // RK stages

    //--------------------------------------------
    // QC
    //--------------------------------------------

    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (verbose>10) fprintf(stdout,"-> check value nan\n");
        //wav_check_value(w_end);
    }

    //--------------------------------------------
    // apply ablexp
    //--------------------------------------------
    if (bdry->is_enable_ablexp == 1) {
       bdry_ablexp_apply(bdry, w_end, wav->ncmp, wav->siz_slice);
    }

    //--------------------------------------------
    // save results
    //--------------------------------------------

    //-- recv by interp
    io_recv_keep(iorecv, w_end, it, wav->ncmp, wav->siz_slice);

    //-- line values
    io_line_keep(ioline, w_end, it, wav->ncmp, wav->siz_slice);

    // snapshot
    if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
      io_snap_nc_put_ac(iosnap, &iosnap_nc, gd, wav, 
                     w_end, w_rhs, nt_total, it, t_end, 1,1,1);
    } else {
      io_snap_nc_put(iosnap, &iosnap_nc, gd, wav, 
                     w_end, w_rhs, nt_total, it, t_end, 1,1,1);
    }

    // zero temp used w_rsh
    wav_zero_edge(gd, wav, w_rhs);

    // swap w_pre and w_end, avoid copying
    w_cur = w_pre; w_pre = w_end; w_end = w_cur;

    if(bdry->is_enable_pml == 1)
    {
      for (int idim=0; idim<CONST_NDIM; idim++) {
        for (int iside=0; iside<2; iside++) {
          bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);
          auxvar->cur = auxvar->pre;
          auxvar->pre = auxvar->end;
          auxvar->end = auxvar->cur;
        }
      }
    }

  } // time loop

  // postproc

  // close nc
  io_snap_nc_close(&iosnap_nc);

  return 0;
}

