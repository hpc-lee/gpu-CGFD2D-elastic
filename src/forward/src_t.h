#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include "constants.h"
#include "gd_t.h"
#include "interp.h"

// cal force_vec_stf/moment_ten_rate 1d index for icmp,it,istage
//  with respect to the start pointer of this source point
#define M_SRC_IND(icmp,it,istage,nt,num_stage) \
  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int total_number;
  int max_nt; // max nt of stf and mrf per src
  int max_stage; // max number of rk stages
  int max_ext; // max extened points

  // for getting value in calculation
  int it;
  int istage;

  // for output
  char evtnm[CONST_MAX_STRLEN];

  // time independent
  int *si; // local i index 
  int *sk; // local k index 
  int *it_begin; // start t index
  int *it_end;   // end   t index
  int   *ext_num; // valid extend points for this src
  int   *ext_indx; // max_ext * total_number
  float *ext_coef;

  // force and/or moment
  int force_actived;
  int moment_actived;

  // time dependent
  // force stf
  float *Fx; // max_stage * max_nt * total_number;
  float *Fz;
  // moment rate
  float *Mxx; // max_stage *max_nt * total_number;
  float *Mzz;
  float *Mxz;
} src_t;

/*************************************************
 * function prototype
 *************************************************/

int
src_read_locate_file(gd_t     *gd,
                     src_t    *src,
                     char     *in_src_file,
                     float     t0,
                     float     dt,
                     int       max_stage,
                     float    *rk_stage_time,
                     int       npoint_half_ext,
                     int       verbose);

float
src_cal_wavelet(float t, char *wavelet_name, float *wavelet_coefs);

float 
fun_ricker(float t, float fc, float t0);

float 
fun_ricker_deriv(float t, float fc, float t0);

float
fun_gauss(float t, float a, float t0);

float
fun_gauss_deriv(float t, float a, float t0);

int
src_set_time(src_t *src, int it, int istage);

int
src_cal_norm_delt2d(float *delt, float x0, float z0,
                    float rx0, float rz0, int LenDelt);

#endif
