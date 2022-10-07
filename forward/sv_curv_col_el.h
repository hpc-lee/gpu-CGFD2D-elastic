#ifndef SV_COL_EL_H
#define SV_COL_EL_H

#include "gd_t.h"
#include "md_t.h"
#include "src_t.h"

/*************************************************
 * function prototype
 *************************************************/

int
sv_curv_col_el_rhs_timg_z2(
             float *restrict  Txx, float *restrict  Tzz,
             float *restrict  Txz, 
             float *restrict hVx , float *restrict hVz ,
             float *restrict xi_x, float *restrict xi_z,
             float *restrict zt_x, float *restrict zt_z,
             float *restrict jac3d, float *restrict slw3d,
             int ni1, int ni2, int nk1, int nk2,
             size_t siz_line,
             int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
             int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
             const int verbose);

int
sv_curv_col_el_rhs_src(
             float *restrict hVx , float *restrict hVz ,
             float *restrict hTxx, float *restrict hTzz,
             float *restrict hTxz, 
             float *restrict jac3d, float *restrict slw3d,
             src_t *src, // short nation for reference member
             const int verbose);

int
sv_curv_graves_Qs(float *w, int ncmp, float dt, gd_t *gd, md_t *md);

#endif
