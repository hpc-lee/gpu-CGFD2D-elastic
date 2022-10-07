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
             float *__restrict__  Txx, float *__restrict__  Tzz,
             float *__restrict__  Txz, 
             float *__restrict__ hVx , float *__restrict__ hVz ,
             float *__restrict__ xi_x, float *__restrict__ xi_z,
             float *__restrict__ zt_x, float *__restrict__ zt_z,
             float *__restrict__ jac3d, float *__restrict__ slw3d,
             int ni1, int ni2, int nk1, int nk2,
             size_t siz_line,
             int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
             int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
             const int verbose);

int
sv_curv_col_el_rhs_src(
             float *__restrict__ hVx , float *__restrict__ hVz ,
             float *__restrict__ hTxx, float *__restrict__ hTzz,
             float *__restrict__ hTxz, 
             float *__restrict__ jac3d, float *__restrict__ slw3d,
             src_t *src, // short nation for reference member
             const int verbose);

int
sv_curv_graves_Qs(float *w, int ncmp, float dt, gd_t *gd, md_t *md);

#endif
