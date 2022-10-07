/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "wav_t.h"

int 
wav_init(gd_t *gd,
         wav_t *V,
         int number_of_levels)
{
  int ierr = 0;

  // Vx,Vz,Txx,Tzz,Txz
  V->ncmp = 5;

  V->nx   = gd->nx;
  V->nz   = gd->nz;
  V->nlevel = number_of_levels;

  V->siz_line   = V->nx;
  V->siz_slice  = V->nx * V->nz;
  V->siz_ilevel = V->siz_slice * V->ncmp;

  // vars
  // 2 Vi, 3 Tij, 4 rk stages
  V->v4d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v4d, wf_el3d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w3d_pos, wf_el3d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w3d_name, wf_el3d_1st");
  
  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_slice;
  }

  // set values
  int icmp = 0;

  /*
   * 0-3: Vx,Vy,Vz
   * 4-9: Txx,Tyy,Tzz,Txz,Tyz,Txy
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txx");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 2;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tzz");
  V->Tzz_pos = cmp_pos[icmp];
  V->Tzz_seq = 3;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txz");
  V->Txz_pos = cmp_pos[icmp];
  V->Txz_seq = 4;
  icmp += 1;

  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int 
wav_ac_init(gd_t *gd,
            wav_t *V,
            int number_of_levels)
{
  int ierr = 0;

  // Vx,Vz,P
  V->ncmp = 3;

  V->nx   = gd->nx;
  V->nz   = gd->nz;
  V->nlevel = number_of_levels;

  V->siz_line   = V->nx;
  V->siz_slice  = V->nx * V->nz;
  V->siz_ilevel = V->siz_slice * V->ncmp;

  // vars
  // 2 Vi, 1 P, 4 rk stages
  V->v4d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v5d, wf_ac3d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w3d_pos, wf_ac3d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w3d_name, wf_ac3d_1st");
  
  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_slice;
  }

  // set values
  int icmp = 0;

  /*
   * 0-1: Vx,Vz
   * 2: P
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","P");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 2;
  icmp += 1;

  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int
wav_check_value(float *__restrict__ w, wav_t *wav)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *ptr = w + icmp * wav->siz_slice;
    for (size_t iptr=0; iptr < wav->siz_slice; iptr++)
    {
      if (ptr[iptr] != ptr[iptr])
      {
        fprintf(stderr, "ERROR: NaN occurs at iptr=%d icmp=%d\n", iptr, icmp);
        fflush(stderr);
        exit(-1);
      }
    }
  }

  return ierr;
}

int
wav_zero_edge(gd_t *gd, wav_t *wav, float *__restrict__ w4d)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *__restrict__ var = w4d + wav->cmp_pos[icmp];

    // z1
    for (int k=0; k < gd->nk1; k++)
    {
      size_t iptr_k = k * gd->siz_line;
        for (int i=0; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // z2
    for (int k=gd->nk2+1; k < gd->nz; k++)
    {
      size_t iptr_k = k * gd->siz_line;
        for (int i=0; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    }

    // x1
    for (int k = gd->nk1; k <= gd->nk2; k++)
    {
      size_t iptr_k = k * gd->siz_line;
        for (int i=0; i < gd->ni1; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

    // x2
    for (int k = gd->nk1; k <= gd->nk2; k++)
    {
      size_t iptr_k = k * gd->siz_line;
        for (int i = gd->ni2+1; i < gd->nx; i++)
        {
          size_t iptr = iptr_k + i;
          var[iptr] = 0.0; 
        }
    } 

  } // icmp

  return ierr;
}
