#ifndef SV_COL_EL_H
#define SV_COL_EL_H

#include "gd_t.h"
#include "md_t.h"
#include "src_t.h"

/*************************************************
 * function prototype
 *************************************************/

__global__ void
sv_curv_col_el_rhs_timg_z2_gpu(
             float *  Txx, float *  Tzz,
             float *  Txz, 
             float * hVx , float * hVz ,
             float * xi_x, float * xi_z,
             float * zt_x, float * zt_z,
             float * jac3d, float * slw3d,
             int ni1, int ni, int nk1, int nk2,
             size_t siz_iz,
             int fdx_len, int * fdx_indx, float * fdx_coef,
             int fdz_len, int * fdz_indx, float * fdz_coef,
             const int verbose);

__global__ void
sv_curv_col_el_rhs_src_gpu(
             float * hVx , float * hVz ,
             float * hTxx, float * hTzz,
             float * hTxz, 
             float * jac3d, float * slw3d,
             src_t src, // short nation for reference member
             const int verbose);

int
sv_curv_graves_Qs(float *w, int ncmp, float dt, gd_t *gd, md_t *md);

#endif
