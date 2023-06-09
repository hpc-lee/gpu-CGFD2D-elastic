/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_el_gpu.h"
#include "sv_curv_col_el_aniso_gpu.h"
#include "cuda_common.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

int
sv_curv_col_el_aniso_onestage(
               float * w_cur_d,
               float * w_rhs_d, 
               wav_t  wav_d,
               fd_wav_t fd_wav_d,
               gd_t   gd_d,
               gdcurv_metric_t  metric_d,
               md_t md_d,
               bdry_t bdry_d,
               src_t src_d,
               // include different order/stentil
               int num_of_fdx_op, fd_op_t *fdx_op,
               int num_of_fdz_op, fd_op_t *fdz_op,
               int fdz_max_len, 
               const int verbose)
{
  // local pointer get each vars
  float * Vx    = w_cur_d + wav_d.Vx_pos ;
  float * Vz    = w_cur_d + wav_d.Vz_pos ;
  float * Txx   = w_cur_d + wav_d.Txx_pos;
  float * Tzz   = w_cur_d + wav_d.Tzz_pos;
  float * Txz   = w_cur_d + wav_d.Txz_pos;
  float * hVx   = w_rhs_d + wav_d.Vx_pos ; 
  float * hVz   = w_rhs_d + wav_d.Vz_pos ; 
  float * hTxx  = w_rhs_d + wav_d.Txx_pos; 
  float * hTzz  = w_rhs_d + wav_d.Tzz_pos; 
  float * hTxz  = w_rhs_d + wav_d.Txz_pos; 

  float * xi_x  = metric_d.xi_x;
  float * xi_z  = metric_d.xi_z;
  float * zt_x  = metric_d.zeta_x;
  float * zt_z  = metric_d.zeta_z;
  float * jac3d = metric_d.jac;

  float * c11   = md_d.c11;
  float * c13   = md_d.c13;
  float * c15   = md_d.c15;
  float * c33   = md_d.c33;
  float * c35   = md_d.c35;
  float * c55   = md_d.c55;
  float * slw3d = md_d.rho;

  // grid size
  int ni1 = gd_d.ni1;
  int ni2 = gd_d.ni2;
  int nk1 = gd_d.nk1;
  int nk2 = gd_d.nk2;

  int ni  = gd_d.ni;
  int nk  = gd_d.nk;
  int nx  = gd_d.nx;
  int nz  = gd_d.nz;
  size_t siz_iz   = gd_d.siz_iz;

  float *vecVx2Vz = bdry_d.vecVx2Vz2;

  // local fd op
  int    fdx_len;
  int    *fdx_indx;
  float  *fdx_coef;
  int    fdz_len;
  int    *fdz_indx;
  float  *fdz_coef;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op
  fdx_len  = fdx_op[num_of_fdx_op-1].total_len;
  fdx_indx = fdx_op[num_of_fdx_op-1].indx;
  fdx_coef = fdx_op[num_of_fdx_op-1].coef;

  fdz_len  = fdz_op[num_of_fdz_op-1].total_len;
  fdz_indx = fdz_op[num_of_fdz_op-1].indx;
  fdz_coef = fdz_op[num_of_fdz_op-1].coef;

  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // allocate max_len because fdz may have different lens
  // these array is for low order surface
  float  fdz_coef_all [num_of_fdz_op*fdz_max_len];
  size_t fdz_shift_all[num_of_fdz_op*fdz_max_len];
  int    fdz_len_all[num_of_fdz_op];
  // loop near surface layers
  for (int n=0; n < num_of_fdz_op-1; n++)
  {
    // get pos and len for this point
    fdz_len_all[n]  = fdz_op[n].total_len;
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_op[n].indx;
    float *p_fdz_coef  = fdz_op[n].coef;
    for (int n_fd = 0; n_fd < fdz_len_all[n] ; n_fd++) {
      fdz_shift_all[n_fd + n*fdz_max_len]  = p_fdz_indx[n_fd] * siz_iz;
      fdz_coef_all [n_fd + n*fdz_max_len]  = p_fdz_coef[n_fd];
    }
  }

  int  *lfdz_len_d = fd_wav_d.fdz_len_d;
  float *lfdx_coef_d = fd_wav_d.fdx_coef_d;
  float *lfdz_coef_d = fd_wav_d.fdz_coef_d;
  float *lfdz_coef_all_d = fd_wav_d.fdz_coef_all_d;
  size_t  *lfdx_shift_d = fd_wav_d.fdx_shift_d;
  size_t  *lfdz_shift_d = fd_wav_d.fdz_shift_d;
  size_t  *lfdz_shift_all_d = fd_wav_d.fdz_shift_all_d;
  int  *lfdx_indx_d = fd_wav_d.fdx_indx_d;
  int  *lfdz_indx_d = fd_wav_d.fdz_indx_d;
  int  *lfdz_indx_all_d = fd_wav_d.fdz_indx_all_d;
  //host to device
  CUDACHECK(cudaMemcpy(lfdx_coef_d,lfdx_coef,fdx_len*sizeof(float),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_coef_d,lfdz_coef,fdz_len*sizeof(float),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdx_shift_d,lfdx_shift,fdx_len*sizeof(size_t),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_shift_d,lfdz_shift,fdz_len*sizeof(size_t),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdx_indx_d,fdx_indx,fdx_len*sizeof(int),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_indx_d,fdz_indx,fdz_len*sizeof(int),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_len_d,fdz_len_all,num_of_fdz_op*sizeof(int),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_coef_all_d,fdz_coef_all,fdz_max_len*num_of_fdz_op*sizeof(float),cudaMemcpyHostToDevice));
  CUDACHECK(cudaMemcpy(lfdz_shift_all_d,fdz_shift_all,fdz_max_len*num_of_fdz_op*sizeof(size_t),cudaMemcpyHostToDevice));

  // inner points
  {
    dim3 block(16,16);
    dim3 grid;
    grid.x = (ni+block.x-1)/block.x;
    grid.y = (nk+block.y-1)/block.y;
    sv_curv_col_el_aniso_rhs_inner_gpu <<<grid, block>>>(
                   Vx,Vz,Txx,Tzz,Txz,
                   hVx,hVz,hTxx,hTzz,hTxz,
                   xi_x, xi_z, zt_x, zt_z,
                   c11,    c13,    c15,    
                           c33,    c35,    
                                   c55,    
                                            slw3d,
                   ni1,ni,nk1,nk,siz_iz,
                   fdx_len, lfdx_shift_d, lfdx_coef_d,
                   fdz_len, lfdz_shift_d, lfdz_coef_d,
                   verbose);
    CUDACHECK( cudaDeviceSynchronize() );
  }

  // free, abs, source in turn
  // free surface at z2
  if (bdry_d.is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // tractiong
    {
      dim3 block(64);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_el_rhs_timg_z2_gpu <<<grid, block>>>(
                     Txx,Tzz,Txz,hVx,hVz,
                     xi_x, xi_z, zt_x, zt_z,
                     jac3d, slw3d,
                     ni1,ni,nk1,nk2,siz_iz,
                     fdx_len, lfdx_indx_d, lfdx_coef_d,
                     fdz_len, lfdz_indx_d, lfdz_coef_d,
                     verbose);
      cudaDeviceSynchronize();
    }

    // velocity: vlow
    {
      dim3 block(64);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_el_aniso_rhs_vlow_z2_gpu <<<grid, block>>> (
                     Vx,Vz,hTxx,hTzz,hTxz,
                     xi_x, xi_z, zt_x, zt_z,
                     c11,    c13,    c15,    
                             c33,    c35,    
                                     c55,    
                                              slw3d,
                     vecVx2Vz,
                     ni1,ni,nk1,nk2,siz_iz,
                     fdx_len, lfdx_shift_d, lfdx_coef_d,
                     num_of_fdz_op,fdz_max_len,lfdz_len_d,
                     lfdz_coef_all_d,lfdz_shift_all_d,
                     verbose);
      cudaDeviceSynchronize();
    }
  }

  // cfs-pml, loop face inside
  if (bdry_d.is_enable_pml == 1)
  {
    sv_curv_col_el_aniso_rhs_cfspml(Vx,Vz,Txx,Tzz,Txz,
                                    hVx,hVz,hTxx,hTzz,hTxz,
                                    xi_x, xi_z, zt_x, zt_z,
                                    c11,    c13,    c15,    
                                            c33,    c35,    
                                                    c55,    
                                                             slw3d,
                                    nk2, siz_iz,
                                    fdx_len, lfdx_shift_d, lfdx_coef_d,
                                    fdz_len, lfdz_shift_d, lfdz_coef_d,
                                    bdry_d,
                                    verbose);
    
  }

  // add source term
  if (src_d.total_number > 0)
  {
    {
      dim3 block(256);
      dim3 grid;
      grid.x = (src_d.total_number+block.x-1) / block.x;
      sv_curv_col_el_rhs_src_gpu  <<< grid,block >>> (
                     hVx,hVz,hTxx,hTzz,hTxz,
                     jac3d, slw3d, 
                     src_d,
                     verbose);
      CUDACHECK( cudaDeviceSynchronize() );
    }
  }

  return 0;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

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
                const int verbose)
{
  // local var
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float slw;
  float c11,    c13,    c15    ;
  float         c33,    c35    ;
  float                 c55    ;
  float xix,xiz,ztx,ztz;

  float * Vx_ptr;
  float * Vz_ptr;
  float * Txx_ptr;
  float * Txz_ptr;
  float * Tzz_ptr;

  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  size_t iz = blockIdx.y * blockDim.y + threadIdx.y;

  // caclu all points
  if(ix<ni && iz<nk)
  {
    size_t iptr = (ix+ni1) + (iz+nk1) * siz_iz;
    Vx_ptr = Vx + iptr;
    Vz_ptr = Vz + iptr;
    Txx_ptr = Txx + iptr;
    Tzz_ptr = Tzz + iptr;
    Txz_ptr = Txz + iptr;

    // Vx derivatives
    M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // Vz derivatives
    M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // Txx derivatives
    M_FD_SHIFT_PTR_MACDRP(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzTxx, Txx_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // Tzz derivatives
    M_FD_SHIFT_PTR_MACDRP(DxTzz, Tzz_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // Txz derivatives
    M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // metric
    xix = xi_x[iptr];
    xiz = xi_z[iptr];
    ztx = zt_x[iptr];
    ztz = zt_z[iptr];

    // medium
    slw = slw3d[iptr];
    c11 = c11d[iptr];
    c13 = c13d[iptr];
    c15 = c15d[iptr];
    c33 = c33d[iptr];
    c35 = c35d[iptr];
    c55 = c55d[iptr];

    // moment equation
    hVx[iptr] = slw*( xix*DxTxx + xiz*DxTxz  
                     +ztx*DzTxx + ztz*DzTxz );
    hVz[iptr] = slw*( xix*DxTxz + xiz*DxTzz 
                     +ztx*DzTxz + ztz*DzTzz );

    // Hooke's equatoin

	  hTxx[iptr] = (c11*xix + c15*xiz) * DxVx + (c15*xix + c13*xiz) * DxVz
               + (c11*ztx + c15*ztz) * DzVx + (c15*ztx + c13*ztz) * DzVz;
    
    hTzz[iptr] = (c13*xix + c35*xiz) * DxVx + (c35*xix + c33*xiz) * DxVz
               + (c13*ztx + c35*ztz) * DzVx + (c35*ztx + c33*ztz) * DzVz;
  
    hTxz[iptr] = (c15*xix + c55*xiz) * DxVx + (c55*xix + c35*xiz) * DxVz
               + (c15*ztx + c55*ztz) * DzVx + (c55*ztx + c35*ztz) * DzVz;
  }

  return;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement vlow boundary
 */

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
                int ni1, int ni, int nk1, int nk2,
                size_t siz_iz, 
                int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
                int num_of_fdz_op, int fdz_max_len, int * fdz_len,
                float *lfdz_coef_all, size_t *lfdz_shift_all,
                const int verbose)
{
  // local var
  int k;
  int n_fd; // loop var for fd
  int lfdz_len;

  // local var
  float DxVx,DxVz;
  float DzVx,DzVz;
  float slw;
  float c11,    c13,    c15    ;
  float         c33,    c35    ;
  float                 c55    ;
  float xix,xiz,ztx,ztz;

  float lfdz_coef[5] = {0.0};
  size_t   lfdz_shift[5] = {0};
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;

  // loop near surface layers
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;
    // get pos and len for this point
    lfdz_len  = fdz_len[n];
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = lfdz_shift_all[n*fdz_max_len+n_fd];
      lfdz_coef [n_fd]  = lfdz_coef_all [n*fdz_max_len+n_fd];
    }

    if(ix<ni)
    {
      size_t iptr   = (ix+ni1) + k * siz_iz;
      // metric
      xix = xi_x[iptr];
      xiz = xi_z[iptr];
      ztx = zt_x[iptr];
      ztz = zt_z[iptr];

      // medium
      slw = slw3d[iptr];
      c11 = c11d[iptr];
      c13 = c13d[iptr];
      c15 = c15d[iptr];
      c33 = c33d[iptr];
      c35 = c35d[iptr];
      c55 = c55d[iptr];

      // Vx derivatives
      M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      // Vz derivatives
      M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      if (k==nk2) // at surface, convert
      {
        size_t ij = (ix+ni1)*4;
        DzVx = vecVx2Vz[ij+2*0+0] * DxVx
             + vecVx2Vz[ij+2*0+1] * DxVz;

        DzVz = vecVx2Vz[ij+2*1+0] * DxVx
             + vecVx2Vz[ij+2*1+1] * DxVz;
      }
      else // lower than surface, lower order
      {
        M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
        M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
      }

      // Hooke's equatoin

	    hTxx[iptr] = (c11*xix + c15*xiz) * DxVx + (c15*xix + c13*xiz) * DxVz
                 + (c11*ztx + c15*ztz) * DzVx + (c15*ztx + c13*ztz) * DzVz;
      
      hTzz[iptr] = (c13*xix + c35*xiz) * DxVx + (c35*xix + c33*xiz) * DxVz
                 + (c13*ztx + c35*ztz) * DzVx + (c35*ztx + c33*ztz) * DzVz;
  
      hTxz[iptr] = (c15*xix + c55*xiz) * DxVx + (c55*xix + c35*xiz) * DxVz
                 + (c15*ztx + c55*ztz) * DzVx + (c55*ztx + c35*ztz) * DzVz;
      
    }
  }

  return;
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

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
               const int verbose)
{
  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip to next face if not cfspml
      if (bdry_d.is_sides_pml[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdry_d.ni1[idim][iside];
      int abs_ni2 = bdry_d.ni2[idim][iside];
      int abs_nk1 = bdry_d.nk1[idim][iside];
      int abs_nk2 = bdry_d.nk2[idim][iside];

      
      int abs_ni = abs_ni2-abs_ni1+1; 
      int abs_nk = abs_nk2-abs_nk1+1; 
      {
        dim3 block(8,8);
        dim3 grid;
        grid.x = (abs_ni+block.x-1)/block.x;
        grid.y = (abs_nk+block.y-1)/block.y;

        sv_curv_col_el_iso_rhs_cfspml_gpu <<<grid, block>>> (
                                idim, iside,
                                Vx, Vz, Txx, Tzz, Txz, 
                                hVx, hVz, hTxx, hTzz, hTxz, 
                                xi_x, xi_z, zt_x, zt_z, 
                                c11d, c13d, c15d, c33d, 
                                c35d, c55d, slw3d,
                                nk2, siz_iz,
                                fdx_len, lfdx_shift, lfdx_coef,
                                fdz_len, lfdz_shift, lfdz_coef,
                                bdry_d, verbose);
        cudaDeviceSynchronize();
      }
    } // iside
  } // idim

  return 0;
}

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
                                  bdry_t bdry_d, const int verbose)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  size_t iz = blockIdx.y * blockDim.y + threadIdx.y;
  float *vecVx2Vz = bdry_d.vecVx2Vz2;

  // local
  size_t iptr, iptr_a;
  float coef_A, coef_B, coef_D, coef_B_minus_1;
  // loop var for fd
  int n_fd; // loop var for fd

  // get index into local var
  int abs_ni1 = bdry_d.ni1[idim][iside];
  int abs_ni2 = bdry_d.ni2[idim][iside];
  int abs_nk1 = bdry_d.nk1[idim][iside];
  int abs_nk2 = bdry_d.nk2[idim][iside];


  int abs_ni = abs_ni2-abs_ni1+1; 
  int abs_nk = abs_nk2-abs_nk1+1; 

  // val on point
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float slw;
  float c11,    c13,    c15    ;
  float         c33,    c35    ;
  float                 c55    ;
  float xix,xiz,ztx,ztz;
  float hVx_rhs,hVz_rhs;
  float hTxx_rhs,hTzz_rhs,hTxz_rhs;
  // for free surface
  float Dx_DzVx,Dx_DzVz;

  // get coef for this face
  float * ptr_coef_A = bdry_d.A[idim][iside];
  float * ptr_coef_B = bdry_d.B[idim][iside];
  float * ptr_coef_D = bdry_d.D[idim][iside];

  bdrypml_auxvar_t *auxvar = &(bdry_d.auxvar[idim][iside]);

  // get pml vars
  float * abs_vars_cur = auxvar->cur;
  float * abs_vars_rhs = auxvar->rhs;

  float * pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
  float * pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
  float * pml_Txx  = abs_vars_cur + auxvar->Txx_pos;
  float * pml_Tzz  = abs_vars_cur + auxvar->Tzz_pos;
  float * pml_Txz  = abs_vars_cur + auxvar->Txz_pos;

  float * pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
  float * pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
  float * pml_hTxx = abs_vars_rhs + auxvar->Txx_pos;
  float * pml_hTzz = abs_vars_rhs + auxvar->Tzz_pos;
  float * pml_hTxz = abs_vars_rhs + auxvar->Txz_pos;

  // for each dim
  if (idim == 0 ) // x direction
  {
    if(ix<abs_ni && iz<abs_nk)
    {
      iptr_a = iz*(abs_ni) + ix;
      iptr   = (ix + abs_ni1) + (iz+abs_nk1) * siz_iz;
      // pml coefs
      // int abs_i = ix;
      coef_D = ptr_coef_D[ix];
      coef_A = ptr_coef_A[ix];
      coef_B = ptr_coef_B[ix];
      coef_B_minus_1 = coef_B - 1.0;

      // metric
      xix = xi_x[iptr];
      xiz = xi_z[iptr];

      // medium
      slw = slw3d[iptr];
      c11 = c11d[iptr];
      c13 = c13d[iptr];
      c15 = c15d[iptr];
      c33 = c33d[iptr];
      c35 = c35d[iptr];
      c55 = c55d[iptr];

      // xi derivatives
      M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      // combine for corr and aux vars
       hVx_rhs = slw * ( xix*DxTxx + xiz*DxTxz );
       hVz_rhs = slw * ( xix*DxTxz + xiz*DxTzz );
      hTxx_rhs = (c11*xix+c15*xiz)*DxVx + (c15*xix+c13*xiz)*DxVz; 
      hTzz_rhs = (c13*xix+c35*xiz)*DxVx + (c35*xix+c33*xiz)*DxVz;
      hTxz_rhs = (c15*xix+c55*xiz)*DxVx + (c55*xix+c35*xiz)*DxVz;

      // 1: make corr to moment equation
      hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
      hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

      // make corr to Hooke's equatoin
      hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
      hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
      hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
      
      // 2: aux var
      //   a1 = alpha + d / beta, dealt in abs_set_cfspml
      pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
      pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
      pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
      pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
      pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];

      // add contributions from free surface condition
      //  not consider timg because conflict with main cfspml,
      //     need to revise in the future if required
      if (bdry_d.is_sides_free[CONST_NDIM-1][1]==1 && (iz+abs_nk1)==nk2)
      {
        // zeta derivatives
        int ij = (ix+abs_ni1)*4;
        Dx_DzVx = vecVx2Vz[ij+2*0+0] * DxVx
                + vecVx2Vz[ij+2*0+1] * DxVz;

        Dx_DzVz = vecVx2Vz[ij+2*1+0] * DxVx
                + vecVx2Vz[ij+2*1+1] * DxVz;

        // metric
        ztx = zt_x[iptr];
        ztz = zt_z[iptr];

        // keep xi derivative terms, including free surface convered
        hTxx_rhs = (c11*ztx+c15*ztz)*Dx_DzVx + (c15*ztx+c13*ztz)*Dx_DzVz; 
        hTzz_rhs = (c13*ztx+c35*ztz)*Dx_DzVx + (c35*ztx+c33*ztz)*Dx_DzVz;
        hTxz_rhs = (c15*ztx+c55*ztz)*Dx_DzVx + (c55*ztx+c35*ztz)*Dx_DzVz;

        // make corr to Hooke's equatoin
        hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
        hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
        hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;

        // aux var
        //   a1 = alpha + d / beta, dealt in abs_set_cfspml
        pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
        pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
        pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
      } // if nk2
    }
  }
  else // z direction
  {
    if(ix<abs_ni && iz<abs_nk)
    {
      iptr_a = iz*abs_ni + ix;
      iptr   = (ix + abs_ni1) + (iz+abs_nk1) * siz_iz;
      // pml coefs
      // int abs_k = iz;
      coef_D = ptr_coef_D[iz];
      coef_A = ptr_coef_A[iz];
      coef_B = ptr_coef_B[iz];
      coef_B_minus_1 = coef_B - 1.0;

      // metric
      ztx = zt_x[iptr];
      ztz = zt_z[iptr];

      // medium
      slw = slw3d[iptr];
      c11 = c11d[iptr];
      c13 = c13d[iptr];
      c15 = c15d[iptr];
      c33 = c33d[iptr];
      c35 = c35d[iptr];
      c55 = c55d[iptr];

      // zt derivatives
      M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

      // combine for corr and aux vars
       hVx_rhs = slw * ( ztx*DzTxx + ztz*DzTxz );
       hVz_rhs = slw * ( ztx*DzTxz + ztz*DzTzz );
      hTxx_rhs = (c11*ztx+c15*ztz)*DzVx + (c15*ztx+c13*ztz)*DzVz; 
      hTzz_rhs = (c13*ztx+c35*ztz)*DzVx + (c35*ztx+c33*ztz)*DzVz;
      hTxz_rhs = (c15*ztx+c55*ztz)*DzVx + (c55*ztx+c35*ztz)*DzVz;

      // 1: make corr to moment equation
      hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
      hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

      // make corr to Hooke's equatoin
      hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
      hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
      hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
      
      // 2: aux var
      //   a1 = alpha + d / beta, dealt in abs_set_cfspml
      pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
      pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
      pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
      pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
      pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];

    } 
  } // if which dim

  return;
}

/*******************************************************************************
 * free surface coef
 * converted matrix for velocity gradient
 *  only implement z2 (top) right now
 ******************************************************************************/

__global__ void
sv_curv_col_el_aniso_dvh2dvz_gpu(gd_t        gd_d,
                                 gdcurv_metric_t metric_d,
                                 md_t       md_d,
                                 bdry_t     bdry_d,
                                 const int verbose)
{
  int ni1 = gd_d.ni1;
  int ni2 = gd_d.ni2;
  int nk1 = gd_d.nk1;
  int nk2 = gd_d.nk2;
  int nx  = gd_d.nx;
  int nz  = gd_d.nz;
  size_t siz_iz   = gd_d.siz_iz;

  // point to each var
  float * xi_x = metric_d.xi_x;
  float * xi_z = metric_d.xi_z;
  float * zt_x = metric_d.zeta_x;
  float * zt_z = metric_d.zeta_z;

  float * c11d = md_d.c11;
  float * c13d = md_d.c13;
  float * c15d = md_d.c15;
  float * c33d = md_d.c33;
  float * c35d = md_d.c35;
  float * c55d = md_d.c55;

  float *vecVx2Vz = bdry_d.vecVx2Vz2;

  float A[2][2], B[2][2], C[2][2];
  float AB[2][2], AC[2][2];

  float c11,    c13,    c15    ;
  float         c33,    c35    ;
  float                 c55    ;
  float xix,xiz,ztx,ztz;
 
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  if(ix<(ni2-ni1+1))
  {
    size_t iptr = (ix+ni1) + nk2 * siz_iz;
    xix = xi_x[iptr];
    xiz = xi_z[iptr];
    ztx = zt_x[iptr];
    ztz = zt_z[iptr];
    
    c11 = c11d[iptr];
    c13 = c13d[iptr];
    c15 = c15d[iptr];
    c33 = c33d[iptr];
    c35 = c35d[iptr];
    c55 = c55d[iptr];

    // first dim: irow; sec dim: jcol, as Fortran code
    A[0][0] = (c11*ztx+c15*ztz)*ztx + (c15*ztx+c55*ztz)*ztz;
    A[1][0] = (c15*ztx+c55*ztz)*ztx + (c13*ztx+c35*ztz)*ztz;

    A[0][1] = (c15*ztx+c13*ztz)*ztx + (c55*ztx+c35*ztz)*ztz; 
    A[1][1] = (c55*ztx+c35*ztz)*ztx + (c35*ztx+c33*ztz)*ztz; 

    fdlib_math_invert2x2(A);
                                                     
    B[0][0] = (c11*xix+c15*xiz)*ztx + (c15*xix+c55*xiz)*ztz;
    B[1][0] = (c15*xix+c55*xiz)*ztx + (c13*xix+c35*xiz)*ztz;

    B[0][1] = (c15*xix+c13*xiz)*ztx + (c55*xix+c35*xiz)*ztz; 
    B[1][1] = (c55*xix+c35*xiz)*ztx + (c35*xix+c33*xiz)*ztz; 
     
    fdlib_math_matmul2x2(A, B, AB);

    size_t ij = (ix+ni1) * 4;

    // save into mat
    for(int irow = 0; irow < 2; irow++){
      for(int jcol = 0; jcol < 2; jcol++){
        vecVx2Vz[ij + irow*2 + jcol] = -AB[irow][jcol];
      }
    }
  }
  return;
}

