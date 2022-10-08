#ifndef WF_EL_1ST_H
#define WF_EL_1ST_H

#include "gd_t.h"

/*************************************************
 * structure
 *************************************************/

/*
 * wavefield structure
 */

// wavefield variables elastic 1st eqn: vel + stress
typedef struct {
  float *v4d; // allocated var

  int n1, n2, n3, n4, n5;
  int nx, nz, ncmp, nlevel;

  size_t siz_iz;
  size_t siz_icmp;
  size_t siz_ilevel;

  size_t *cmp_pos;
  char  **cmp_name;
  size_t *level_pos;

  size_t Vx_pos;
  size_t Vz_pos;
  size_t Txx_pos;
  size_t Tzz_pos;
  size_t Txz_pos;

  // sequential index 0-based
  size_t Vx_seq;
  size_t Vz_seq;
  size_t Txx_seq;
  size_t Tzz_seq;
  size_t Txz_seq;

} wav_t;

/*************************************************
 * function prototype
 *************************************************/

int 
wav_init(gd_t *gd,
         wav_t *V,
         int number_of_levels);

int
wav_check_value(float *__restrict__ w, wav_t *wav);

int
wav_zero_edge(gd_t *gd, wav_t *wav,
              float *__restrict__ w4d);

int 
wav_ac_init(gd_t *gd,
            wav_t *V,
            int number_of_levels);

#endif
