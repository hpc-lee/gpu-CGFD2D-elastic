#ifndef BDRY_H
#define BDRY_H

#include "constants.h"
#include "gd_t.h"
#include "wav_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  float *var;
  int nx, nz, ncmp, nlevel;

  size_t siz_line;
  size_t siz_slice;
  size_t siz_ilevel;

  size_t *cmp_pos;
  char  **cmp_name;
  size_t *level_pos;

  size_t Vx_pos;
  size_t Vz_pos;
  size_t Txx_pos;
  size_t Tzz_pos;
  size_t Txz_pos;

  // for rk scheme
  float *cur;
  float *pre;
  float *rhs;
  float *end;
  float *tmp;

} bdrypml_auxvar_t;

/*
 * structure for block index range
 */

typedef struct {
  int enable;
  int nx1,nx2,nz1,nz2,nx,nz;
  int ni1,ni2,nk1,nk2,ni,nk;
} bdry_block_t;

typedef struct
{
  // 0 or 1 to indicate corresponding boundary condition
  //   such implementation simplifies calling the funcs
  int is_sides_pml [CONST_NDIM][2];
  int is_sides_free[CONST_NDIM][2];
  int is_sides_mpml[CONST_NDIM][2];
  int is_sides_ablexp [CONST_NDIM][2];

  int is_enable_pml;
  int is_enable_free;
  int is_enable_mpml;
  int is_enable_ablexp;

  // same as grid, to make here self contained
  int nx;
  int nz;

  // used for PML or exp
  int num_of_layers[CONST_NDIM][2]; //

  int ni1[CONST_NDIM][2];
  int ni2[CONST_NDIM][2];
  int nk1[CONST_NDIM][2];
  int nk2[CONST_NDIM][2];

  //
  // for ADE CFS-PML
  //

  float *A[CONST_NDIM][2]; // dim, side, length
  float *B[CONST_NDIM][2]; // dim, side, length
  float *D[CONST_NDIM][2]; // dim, side, length

  bdrypml_auxvar_t auxvar[CONST_NDIM][2];

  //
  // for ABLEXP
  //

  // use 4 blocks to partition the boundaries
  bdry_block_t bdry_blk[CONST_NDIM_2];

  float *ablexp_Ex;
  float *ablexp_Ez;

  // top
  float *vecVx2Vz2; // [j,i, dzVi, dxVi]

  // bottom
  float *vecVx2Vz1;

  // left
  float *vecVy2Vx1;

  // right
  float *vecVy2Vx2;

  // front
  float *vecVx2Vy1;

  // back
  float *vecVx2Vy2;

} bdry_t;

/*************************************************
 * function prototype
 *************************************************/

int
bdry_init(bdry_t *bdry, int nx, int nz);

int
bdry_free_set(gd_t    *gd,
              bdry_t  *bdryfree,
              int   in_is_sides[][2],
              const int verbose);

int
bdry_pml_set(gd_t     *gd,
             wav_t    *wav,
             bdry_t   *bdrypml,
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_alpha_max[][2], 
             float in_beta_max[][2], 
             float in_velocity[][2], 
             int verbose);

int
bdry_pml_auxvar_init(int nx, int nz, 
                     wav_t *wav,
                     bdrypml_auxvar_t *auxvar,
                     const int verbose);

float
bdry_pml_cal_R(float N);

float
bdry_pml_cal_dmax(float L, float Vp, float Rpp);

float
bdry_pml_cal_amax(float fc);

float
bdry_pml_cal_d(float x, float L, float dmax);

float
bdry_pml_cal_a(float x, float L, float amax);

float
bdry_pml_cal_b(float x, float L, float bmax);

int
bdry_cal_abl_len_dh(gd_t *gd, 
                    int abs_ni1, int abs_ni2,
                    int abs_nk1, int abs_nk2,
                    int idim,
                    float *avg_L, float *avg_dh);

int
bdry_ablexp_set(gd_t *gd,
                wav_t *wav,
                bdry_t *bdryexp,
                int   in_is_sides[][2],
                int   in_num_layers[][2],
                float in_velocity[][2], //
                float dt,
                int verbose);

float
bdry_ablexp_cal_mask(int i, float vel, float dt, int num_lay, float dh);

int
bdry_ablexp_apply(bdry_t *bdry, float *w_end, int ncmp, size_t siz_slice);

#endif
