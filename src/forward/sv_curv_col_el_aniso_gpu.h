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
               float * w_cur_d,
               float * w_rhs_d, 
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
sv_curv_col_el_aniso_rhs_inner_gpu(
                float *  Vx , float *  Vz ,
                float *  Txx, float *  Tzz,
                float *  Txz, 
                float * hVx , float * hVz ,
                float * hTxx, float * hTzz,
                float * hTxz, 
                float * xi_x, float * xi_z,
                float * zt_x, float * zt_z,
                float * c11d, float * c13d,
                float * c15d, float * c33d,
                float * c35d, float * c55d,
                float * slw3d,
                int ni1, int ni, int nk1, int nk,
                size_t siz_iz,
                int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
                int fdz_len, size_t * lfdz_shift, float * lfdz_coef,
                const int verbose);

__global__ void
sv_curv_col_el_aniso_rhs_vlow_z2_gpu(
                float *  Vx , float *  Vz ,
                float * hTxx, float * hTzz,
                float * hTxz, 
                float * xi_x, float * xi_z,
                float * zt_x, float * zt_z,
                float * c11d, float * c13d,
                float * c15d, float * c33d,
                float * c35d, float * c55d,
                float * slw3d,
                float * vecVx2Vz,
               int ni1, int ni, int nk1, int nk,
               size_t siz_iz, 
               int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
               int num_of_fdz_op, int fdz_max_len, int * fdz_len,
               float *lfdz_coef_all, size_t *lfdz_shift_all,
                const int verbose);

int
sv_curv_col_el_aniso_rhs_cfspml(
               float *  Vx , float *  Vz ,
               float *  Txx, float *  Tzz,
               float *  Txz, 
               float * hVx , float * hVz ,
               float * hTxx, float * hTzz,
               float * hTxz, 
               float * xi_x, float * xi_z,
               float * zt_x, float * zt_z,
               float * c11d, float * c13d,
               float * c15d, float * c33d,
               float * c35d, float * c55d,
               float * slw3d,
               int nk2, size_t siz_iz,
               int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
               int fdz_len, size_t * lfdz_shift, float * lfdz_coef,
               bdry_t bdry_d,
               const int verbose);

__global__ void
sv_curv_col_el_iso_rhs_cfspml_gpu(int idim, int iside,
                                  float * Vx,    float * Vz,
                                  float * Txx,   float * Tzz, 
                                  float * Txz, 
                                  float * hVx,   float * hVz,
                                  float * hTxx,  float * hTzz, 
                                  float * hTxz, 
                                  float * xi_x,  float * xi_z,
                                  float * zt_x,  float * zt_z,
                                  float * c11d, float * c13d, 
                                  float * c15d, float * c33d, 
                                  float * c35d, float * c55d, 
                                  float * slw3d,
                                  int nk2, size_t siz_iz,
                                  int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
                                  int fdz_len, size_t * lfdz_shift, float * lfdz_coef,
                                  bdry_t bdry_d, const int verbose);

__global__ void
sv_curv_col_el_aniso_dvh2dvz_gpu(gd_t        gd_d,
                                 gdcurv_metric_t metric_d,
                                 md_t       md_d,
                                 bdry_t     bdry_d,
                                 const int verbose);

#endif
