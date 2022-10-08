#ifndef SV_CURV_COL_AC_H
#define SV_CURV_COL_AC_H

#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

int
sv_curv_col_ac_iso_onestage(
            float *__restrict__ w_cur_d,
            float *__restrict__ w_rhs_d, 
            wav_t  wav_d,
            fd_wav_t  fd_wav_d,
            gd_t   gd_d,
            gdcurv_metric_t  metric_d,
            md_t md_d,
            bdry_t bdry_d,
            src_t src_d,
            // include different order/stentil
            int num_of_fdx_op, fd_op_t *fdx_op,
            int num_of_fdz_op, fd_op_t *fdz_op,
            int fdz_max_len, 
            const int verbose);

__global__ void
sv_curv_col_ac_iso_rhs_inner(
              float *__restrict__  Vx , float *__restrict__  Vz ,
              float *__restrict__  P, 
              float *__restrict__ hVx , float *__restrict__ hVz ,
              float *__restrict__ hP, 
              float *__restrict__ xi_x, float *__restrict__ xi_z,
              float *__restrict__ zt_x, float *__restrict__ zt_z,
              float *__restrict__ kappa3d, float *__restrict__ slw3d,
              int ni1, int ni, int nk1, int nk,
              size_t siz_iz,
              int fdx_len, int *__restrict__ lfdx_shift, float *__restrict__ lfdx_coef,
              int fdz_len, int *__restrict__ lfdz_shift, float *__restrict__ lfdz_coef,
              const int verbose);

__global__ void
sv_curv_col_ac_iso_rhs_timg_z2_gpu(
               float *__restrict__  P,
               int ni1, int ni, int nk1, int nk, int nz,
               size_t siz_iz,
               const int verbose);

__global__ void
sv_curv_col_ac_iso_rhs_vlow_z2_gpu(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__ hP, 
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ kappa3d, float *__restrict__ slw3d,
               int ni1, int ni, int nk1, int nk,
               size_t siz_iz, 
               int fdx_len, int *__restrict__ lfdx_shift, float *__restrict__ lfdx_coef,
               int num_of_fdz_op, int fdz_max_len, int * fdz_len,
               float *lfdz_coef_all, size_t *lfdz_shift_all,
               const int verbose);

int
sv_curv_col_ac_iso_rhs_cfspml(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__  P, 
               float *__restrict__ hVx , float *__restrict__ hVz ,
               float *__restrict__ hP,
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ kappa3d, float *__restrict__ slw3d,
               int nk2, size_t siz_iz,
               int fdx_len, int *__restrict__ lfdx_shift, float *__restrict__ lfdx_coef,
               int fdz_len, int *__restrict__ lfdz_shift, float *__restrict__ lfdz_coef,
               bdry_t bdry_d,
               const int verbose);

__global__ void
sv_curv_col_ac_iso_rhs_src_gpu(
             float *__restrict__ hVx , float *__restrict__ hVz ,
             float *__restrict__ hP, 
             float *__restrict__ jac3d, float *__restrict__ slw3d,
             src_t src_d, // short nation for reference member
             const int verbose);

#endif
