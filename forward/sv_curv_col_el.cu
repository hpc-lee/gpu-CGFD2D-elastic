#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"

#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
/*
 * implement traction image boundary 
 */

int
sv_curv_col_el_rhs_timg_z2(
                    float *__restrict__  Txx, float *__restrict__  Tzz,
                    float *__restrict__  Txz, 
                    float *__restrict__ hVx , float *__restrict__ hVz ,
                    float *__restrict__ xi_x, float *__restrict__ xi_z,
                    float *__restrict__ zt_x, float *__restrict__ zt_z,
                    float *__restrict__ jac3d, float *__restrict__ slw3d,
                    int ni1, int ni2, int nk1, int nk2,
                    size_t siz_iz, 
                    int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
                    int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
                    const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTx,DzTz;
  float slwjac;
  float xix,xiz,ztx,ztz;

  // to save traction and other two dir force var
  float vecxi[fdx_len];
  float veczt[fdz_len];
  int n, iptr4vec;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // last indx, free surface force Tx/Ty/Tz to 0 in cal
  int k_min = nk2 - fdz_indx[fdz_len-1];

  // point affected by timg
  for (size_t k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    int n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
          // metric
          xix = xi_x[iptr];
          xiz = xi_z[iptr];
          ztx = zt_x[iptr];
          ztz = zt_z[iptr];

          // slowness and jac
          slwjac = slw3d[iptr] / jac3d[iptr];

          //
          // for hVx
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txx[iptr4vec]
                                          + xi_z[iptr4vec] * Txz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_iz;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                          + zt_z[iptr4vec] * Txz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++)
          {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            //int n_img = n_free-(n-n_free);
            iptr4vec = iptr + n_img * siz_iz;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                           + zt_z[iptr4vec] * Txz[iptr4vec] );
            //veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DzTz ) * slwjac;

          //
          // for hVz
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txz[iptr4vec]
                                          + xi_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_iz;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                          + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_iz;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                           + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // for hVx 
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVz[iptr] = ( DxTx+DzTz ) * slwjac;

          // next
          iptr += 1;
      }
  }

  return 0;
}

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_curv_col_el_rhs_src(
            float *__restrict__ hVx , float *__restrict__ hVz ,
            float *__restrict__ hTxx, float *__restrict__ hTzz,
            float *__restrict__ hTxz, 
            float *__restrict__ jac3d, float *__restrict__ slw3d,
            src_t *src, // short nation for reference member
            const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
  float fx, fz;
  float Mxx,Mzz,Mxz;

  int it     = src->it;
  int istage = src->istage;

  // add src; is is a commont iterater var
  for (int is=0; is < src->total_number; is++)
  {
    int   it_start = src->it_begin[is];
    int   it_end   = src->it_end  [is];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src->ext_indx + is * max_ext;
      float *ptr_ext_coef = src->ext_coef + is * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   is * src->max_nt * src->max_stage // skip other src
                           + it_to_it_start * src->max_stage // skip other time step
                           + istage;
      if (src->force_actived == 1) {
        fx  = src->Fx [iptr_cur_stage];
        fz  = src->Fz [iptr_cur_stage];
      }
      if (src->moment_actived == 1) {
        Mxx = src->Mxx[iptr_cur_stage];
        Mzz = src->Mzz[iptr_cur_stage];
        Mxz = src->Mxz[iptr_cur_stage];
      }
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        //fprintf(stdout,"-----> src_iptr=%d,src_ceof=%g\n",iptr,coef);

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          hVx[iptr] += fx * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          hTxx[iptr] -= Mxx * rjac;
          hTzz[iptr] -= Mzz * rjac;
          hTxz[iptr] -= Mxz * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

int
sv_curv_graves_Qs(float *w, int ncmp, float dt, gd_t *gd, md_t *md)
{
  int ierr = 0;

  float coef = - PI * md->visco_Qs_freq * dt;

  for (int icmp=0; icmp<ncmp; icmp++)
  {
    float *__restrict__ var = w + icmp * gd->siz_icmp;

    for (int k = gd->nk1; k <= gd->nk2; k++)
    {
        for (int i = gd->ni1; i <= gd->ni2; i++)
        {
          size_t iptr = i + k * gd->siz_iz;

          float Qatt = expf( coef / md->Qs[iptr] );

          var[iptr] *= Qatt;
        }
    }
  }

  return ierr;
}
