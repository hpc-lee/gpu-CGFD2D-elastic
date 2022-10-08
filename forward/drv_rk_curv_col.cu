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
#include "sv_curv_col_el_gpu.h"
#include "sv_curv_col_el_iso_gpu.h"
#include "sv_curv_col_el_vti_gpu.h"
#include "sv_curv_col_el_aniso_gpu.h"
#include "sv_curv_col_ac_iso_gpu.h"
#include "alloc.h"
#include "cuda_common.h"

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
  int fdz_max_len      = fd->fdz_max_len;
  int num_of_fdz_op    = fd->num_of_fdz_op;

  //gdinfo
  int ni = gd->ni;
  int nk = gd->nk;

  // local allocated array
  char ou_file[CONST_MAX_STRLEN];

  gd_t   gd_d;
  md_t   md_d;
  src_t  src_d;
  wav_t  wav_d;
  bdry_t  bdry_d;
  gdcurv_metric_t metric_d;
  fd_wav_t fd_wav_d;

  // init device struct, and copy data from host to device
  init_gdinfo_device(gd, &gd_d);
  init_md_device(md, &md_d);
  init_fd_device(fd, &fd_wav_d);
  init_src_device(src, &src_d);
  init_metric_device(metric, &metric_d);
  init_bdry_device(gd, bdry, &bdry_d);
  init_wave_device(wav, &wav_d);
  //---------------------------------------
  // get device wavefield 
  float *__restrict__ w_buff = wav->v4d; // size number is V->siz_icmp * (V->ncmp+6)

  // local pointer
  float *__restrict__ w_cur_d;
  float *__restrict__ w_pre_d;
  float *__restrict__ w_rhs_d;
  float *__restrict__ w_end_d;
  float *__restrict__ w_tmp_d;

  // get wavefield
  w_pre_d = wav_d.v4d + wav_d.siz_ilevel * 0; // previous level at n
  w_tmp_d = wav_d.v4d + wav_d.siz_ilevel * 1; // intermidate value
  w_rhs_d = wav_d.v4d + wav_d.siz_ilevel * 2; // for rhs
  w_end_d = wav_d.v4d + wav_d.siz_ilevel * 3; // end level at n+1

  int   ipair, istage;
  float t_cur;
  float t_end; // time after this loop for nc output

  // create snapshot nc output files
  if (verbose>0) fprintf(stdout,"prepare snap nc output ...\n"); 
  iosnap_nc_t  iosnap_nc;
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    io_snap_nc_create_ac(iosnap, &iosnap_nc);
  } else {
    io_snap_nc_create(iosnap, &iosnap_nc);
  }
  // set pml for rk
  if(bdry_d.is_enable_pml == 1)
  {
    for (int idim=0; idim<CONST_NDIM; idim++) {
      for (int iside=0; iside<2; iside++) {
        if (bdry_d.is_sides_pml[idim][iside]==1) {
          bdrypml_auxvar_t *auxvar_d = &(bdry_d.auxvar[idim][iside]);
          auxvar_d->pre = auxvar_d->var + auxvar_d->siz_ilevel * 0;
          auxvar_d->tmp = auxvar_d->var + auxvar_d->siz_ilevel * 1;
          auxvar_d->rhs = auxvar_d->var + auxvar_d->siz_ilevel * 2;
          auxvar_d->end = auxvar_d->var + auxvar_d->siz_ilevel * 3;
        }
      }
    }
  }

  // calculate conversion matrix for free surface
  if (bdry_d.is_sides_free[CONST_NDIM-1][1] == 1)
  {
    if(md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
    {
      dim3 block(128);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_el_iso_dvh2dvz_gpu(gd_d, metric_d, md_d, bdry_d, verbose);
      CUDACHECK(cudaDeviceSynchronize());
    } 
    else if(md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) 
    {
      dim3 block(128);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_el_aniso_dvh2dvz_gpu(gd_d, metric_d, md_d, bdry_d, verbose);
      CUDACHECK(cudaDeviceSynchronize());
    }
    else if(md->medium_type == CONST_MEDIUM_ELASTIC_VTI)
    {
      dim3 block(128);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_el_vti_dvh2dvz_gpu(gd_d, metric_d, md_d, bdry_d, verbose);
      CUDACHECK(cudaDeviceSynchronize());
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
        w_cur_d = w_pre_d;
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              bdry_d.auxvar[idim][iside].cur = bdry_d.auxvar[idim][iside].pre;
            }
          }
        }
      }
      else
      {
        w_cur_d = w_tmp_d;
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              bdry_d.auxvar[idim][iside].cur = bdry_d.auxvar[idim][iside].tmp;
            }
          }
        }
      }

      // set src_t time
      src_set_time(&src_d, it, istage);

      // compute rhs
      switch (md_d.medium_type)
      {
        case CONST_MEDIUM_ELASTIC_ISO : {
          sv_curv_col_el_iso_onestage(
              w_cur_d,w_rhs_d,wav_d,fd_wav_d,
              gd_d, metric_d, md_d, bdry_d, src_d,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
          break;
        }

        case CONST_MEDIUM_ELASTIC_VTI : {
          sv_curv_col_el_vti_onestage(
              w_cur_d,w_rhs_d,wav_d,fd_wav_d,
              gd_d, metric_d, md_d, bdry_d, src_d,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
          break;
        }

        case CONST_MEDIUM_ELASTIC_ANISO : {
          sv_curv_col_el_aniso_onestage(
              w_cur_d,w_rhs_d,wav_d,fd_wav_d,
              gd_d, metric_d, md_d, bdry_d, src_d,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
          break;
        }

        case CONST_MEDIUM_ACOUSTIC_ISO : {
          sv_curv_col_ac_iso_onestage(
              w_cur_d,w_rhs_d,wav_d,fd_wav_d,
              gd_d, metric_d, md_d, bdry_d, src_d,
              fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
              fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
              fd->fdz_max_len,
              verbose);
          break;
        }
        //  synchronize onestage device func.
        CUDACHECK(cudaDeviceSynchronize());
      }

      // rk start
      if (istage==0)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // wavefield
        {
          dim3 block(256);
          dim3 grid;
          grid.x = (wav_d.siz_ilevel + block.x - 1) / block.x;
          wav_update <<<grid, block>>> (wav_d.siz_ilevel, coef_a, w_tmp_d, w_pre_d, w_rhs_d);
        }

        // apply Qs
        //if (md->visco_type == CONST_VISCO_GRAVES) {
        //  sv_curv_graves_Qs(w_tmp, wave->ncmp, gd, md);
        //}

        // pml_tmp
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry_d.is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
                dim3 block(256);
                dim3 grid;
                grid.x = (auxvar_d->siz_ilevel + block.x - 1) / block.x;
                wav_update <<<grid, block>>> (
                           auxvar_d->siz_ilevel, coef_a, auxvar_d->tmp, auxvar_d->pre, auxvar_d->rhs);
              }
            }
          }
        }

        // w_end
        {
          dim3 block(256);
          dim3 grid;
          grid.x = (wav_d.siz_ilevel + block.x - 1) / block.x;
          wav_update <<<grid, block>>> (wav_d.siz_ilevel, coef_b, w_end_d, w_pre_d, w_rhs_d);
        }
        // pml_end
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry_d.is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
                dim3 block(256);
                dim3 grid;
                grid.x = (auxvar_d->siz_ilevel + block.x - 1) / block.x;
                wav_update <<<grid, block>>> (
                           auxvar_d->siz_ilevel, coef_b, auxvar_d->end, auxvar_d->pre, auxvar_d->rhs);
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
        {
          dim3 block(256);
          dim3 grid;
          grid.x = (wav_d.siz_ilevel + block.x - 1) / block.x;
          wav_update <<<grid, block>>> (wav_d.siz_ilevel, coef_a, w_tmp_d, w_pre_d, w_rhs_d);
        }

        // apply Qs
        //if (md->visco_type == CONST_VISCO_GRAVES) {
        //  sv_curv_graves_Qs(w_tmp, wave->ncmp, gd, md);
        //}

        // pml_tmp
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry_d.is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
                dim3 block(256);
                dim3 grid;
                grid.x = (auxvar_d->siz_ilevel + block.x - 1) / block.x;
                wav_update <<<grid, block>>> (
                           auxvar_d->siz_ilevel, coef_a, auxvar_d->tmp, auxvar_d->pre, auxvar_d->rhs);
              }
            }
          }
        }

        // w_end
        {
          dim3 block(256);
          dim3 grid;
          grid.x = (wav_d.siz_ilevel + block.x - 1) / block.x;
          wav_update_end <<<grid, block>>> (wav_d.siz_ilevel, coef_b, w_end_d, w_rhs_d);
        }
        // pml_end
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry_d.is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
                dim3 block(256);
                dim3 grid;
                grid.x = (auxvar_d->siz_ilevel + block.x - 1) / block.x;
                wav_update_end <<<grid, block>>> (
                           auxvar_d->siz_ilevel, coef_b, auxvar_d->end, auxvar_d->rhs);
              }
            }
          }
        }
      }
      else // last stage
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        {
          dim3 block(256);
          dim3 grid;
          grid.x = (wav_d.siz_ilevel + block.x - 1) / block.x;
          wav_update_end <<<grid, block>>>(wav_d.siz_ilevel, coef_b, w_end_d, w_rhs_d);
        }

        // apply Qs
        //if (md->visco_type == CONST_VISCO_GRAVES_QS) {
        //  sv_curv_graves_Qs(w_end, wav->ncmp, dt, gd, md);
        //}
        
        // pml_end
        if(bdry_d.is_enable_pml == 1)
        {
          for (int idim=0; idim<CONST_NDIM; idim++) {
            for (int iside=0; iside<2; iside++) {
              if (bdry_d.is_sides_pml[idim][iside]==1) {
                bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
                dim3 block(256);
                dim3 grid;
                grid.x = (auxvar_d->siz_ilevel + block.x - 1) / block.x;
                wav_update_end <<<grid, block>>> (
                           auxvar_d->siz_ilevel, coef_b, auxvar_d->end, auxvar_d->rhs);
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
    if (bdry_d.is_enable_ablexp == 1) {
       bdry_ablexp_apply(bdry_d, gd, w_end_d, wav->ncmp);
    }

    //--------------------------------------------
    // save results
    //--------------------------------------------

    //-- recv by interp
    io_recv_keep(iorecv, w_end_d, w_buff, it, wav->ncmp, wav->siz_icmp);

    //-- line values
    io_line_keep(ioline, w_end_d, w_buff, it, wav->ncmp, wav->siz_icmp);

    // snapshot
    if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
      io_snap_nc_put_ac(iosnap, &iosnap_nc, gd, wav, 
                     w_end_d, w_buff, nt_total, it, t_end, 1,1,1);
    } else {
      io_snap_nc_put(iosnap, &iosnap_nc, gd, wav, 
                     w_end_d, w_buff, nt_total, it, t_end, 1,1,1);
    }

    // swap w_pre and w_end, avoid copying
    w_cur_d = w_pre_d; w_pre_d = w_end_d; w_end_d = w_cur_d;

    if(bdry_d.is_enable_pml == 1)
    {
      for (int idim=0; idim<CONST_NDIM; idim++) {
        for (int iside=0; iside<2; iside++) {
          bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);
          auxvar_d->cur = auxvar_d->pre;
          auxvar_d->pre = auxvar_d->end;
          auxvar_d->end = auxvar_d->cur;
        }
      }
    }

  } // time loop

  dealloc_md_device(md_d);
  dealloc_fd_device(fd_wav_d);
  dealloc_metric_device(metric_d);
  dealloc_src_device(src_d);
  dealloc_bdry_device(bdry_d);
  dealloc_wave_device(wav_d);

  // postproc

  // close nc
  io_snap_nc_close(&iosnap_nc);

  return 0;
}

