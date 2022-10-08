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
             float *__restrict__ w_cur,
             float *__restrict__ rhs, 
             wav_t  *wav,
             gd_t   *gd,
             gdcurv_metric_t  *metric,
             md_t *md,
             bdry_t *bdry,
             src_t *src,
             // include different order/stentil
             int num_of_fdx_op, fd_op_t *fdx_op,
             int num_of_fdz_op, fd_op_t *fdz_op,
             int fdz_max_len, 
             const int verbose);

int
sv_curv_col_ac_iso_rhs_inner(
             float *__restrict__  Vx , float *__restrict__  Vz ,
             float *__restrict__  P, 
             float *__restrict__ hVx , float *__restrict__ hVz ,
             float *__restrict__ hP, 
             float *__restrict__ xi_x, float *__restrict__ xi_z,
             float *__restrict__ zt_x, float *__restrict__ zt_z,
             float *__restrict__ kappa3d, float *__restrict__ slw3d,
             int ni1, int ni2, int nk1, int nk2,
             size_t siz_iz,
             int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
             int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
             const int verbose);

int
sv_curv_col_ac_iso_rhs_timg_z2(
             float *__restrict__  P,
             int ni1, int ni2, int nk1, int nk2, int nz,
             size_t siz_iz, 
             const int verbose);

int
sv_curv_col_ac_iso_rhs_vlow_z2(
             float *__restrict__  Vx , float *__restrict__  Vz ,
             float *__restrict__ hP, 
             float *__restrict__ xi_x, float *__restrict__ xi_z,
             float *__restrict__ zt_x, float *__restrict__ zt_z,
             float *__restrict__ kappa3d, float *__restrict__ slw3d,
             int ni1, int ni2, int nk1, int nk2,
             size_t siz_iz,
             int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
             int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
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
               int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
               int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
               bdry_t *bdry,
               const int verbose);

int
sv_curv_col_ac_iso_rhs_src(
               float *__restrict__ hVx , float *__restrict__ hVz ,
               float *__restrict__ hP, 
               float *__restrict__ jac3d, float *__restrict__ slw3d,
               src_t *src, // short nation for reference member
               const int verbose);

#endif
