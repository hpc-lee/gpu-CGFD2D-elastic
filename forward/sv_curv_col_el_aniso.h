#ifndef SV_CURV_COL_EL_ANISO_H
#define SV_CURV_COL_EL_ANISO_H

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
sv_curv_col_el_aniso_onestage(
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
sv_curv_col_el_aniso_rhs_inner(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__  Txx, float *__restrict__  Tzz,
               float *__restrict__  Txz, 
               float *__restrict__ hVx , float *__restrict__ hVz ,
               float *__restrict__ hTxx, float *__restrict__ hTzz,
               float *__restrict__ hTxz, 
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ c11d, float *__restrict__ c13d,
               float *__restrict__ c15d, float *__restrict__ c33d,
               float *__restrict__ c35d, float *__restrict__ c55d,
               float *__restrict__ slw3d,
               int ni1, int ni2, int nk1, int nk2,
               size_t siz_line,
               int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
               int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
               const int verbose);

int
sv_curv_col_el_aniso_rhs_vlow_z2(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__ hTxx, float *__restrict__ hTzz,
               float *__restrict__ hTxz, 
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ c11d, float *__restrict__ c13d,
               float *__restrict__ c15d, float *__restrict__ c33d,
               float *__restrict__ c35d, float *__restrict__ c55d,
               float *__restrict__ slw3d,
               float *__restrict__ vecVx2Vz,
               int ni1, int ni2, int nk1, int nk2,
               size_t siz_line,
               int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
               int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
               const int verbose);

int
sv_curv_col_el_aniso_rhs_cfspml(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__  Txx, float *__restrict__  Tzz,
               float *__restrict__  Txz, 
               float *__restrict__ hVx , float *__restrict__ hVz ,
               float *__restrict__ hTxx, float *__restrict__ hTzz,
               float *__restrict__ hTxz, 
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ c11d, float *__restrict__ c13d,
               float *__restrict__ c15d, float *__restrict__ c33d,
               float *__restrict__ c35d, float *__restrict__ c55d,
               float *__restrict__ slw3d,
               int nk2, size_t siz_line,
               int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
               int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
               bdry_t *bdry,
               const int verbose);

int
sv_curv_col_el_aniso_dvh2dvz(gd_t        *gd,
                             gdcurv_metric_t *metric,
                             md_t       *md,
                             bdry_t     *bdry,
                             const int verbose);

#endif
