/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_ac_iso_gpu.h"
#include "cuda_common.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

int
sv_curv_col_ac_iso_onestage(
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
            const int verbose)
{
  // local pointer get each vars
  float * Vx    = w_cur_d + wav_d.Vx_pos ;
  float * Vz    = w_cur_d + wav_d.Vz_pos ;
  float * P     = w_cur_d + wav_d.Txx_pos;
  float * hVx   = w_rhs_d + wav_d.Vx_pos ; 
  float * hVz   = w_rhs_d + wav_d.Vz_pos ; 
  float * hP    = w_rhs_d + wav_d.Txx_pos; 

  float * xi_x  = metric_d.xi_x;
  float * xi_z  = metric_d.xi_z;
  float * zt_x  = metric_d.zeta_x;
  float * zt_z  = metric_d.zeta_z;
  float * jac3d = metric_d.jac;

  float * kappa3d = md_d.kappa;
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

  // free surface at z2 for pressure
  if (bdry_d.is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // imaging
    {
      dim3 block(64);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_ac_iso_rhs_timg_z2_gpu <<<grid, block>>>(
                                     P,
                                     ni1,ni,nk1,nk2,nz,
                                     siz_iz,
                                     verbose);
      cudaDeviceSynchronize();
    }
  }

  // inner points
  {
    dim3 block(16,16);
    dim3 grid;
    grid.x = (ni+block.x-1)/block.x;
    grid.y = (nk+block.y-1)/block.y;
    sv_curv_col_ac_iso_rhs_inner_gpu <<<grid, block>>>(
                Vx,Vz,P,
                hVx,hVz,hP,
                xi_x, xi_z, zt_x, zt_z,
                kappa3d, slw3d,
                ni1,ni,nk1,nk,siz_iz,
                fdx_len, lfdx_shift_d, lfdx_coef_d,
                fdz_len, lfdz_shift_d, lfdz_coef_d,
                verbose);
    CUDACHECK( cudaDeviceSynchronize() );
  }

  // free, abs, source in turn
  // free surface at z2 for velocity
  if (bdry_d.is_sides_free[CONST_NDIM-1][1] == 1)
  {
    {
      dim3 block(64);
      dim3 grid;
      grid.x = (ni+block.x-1)/block.x;
      sv_curv_col_ac_iso_rhs_vlow_z2_gpu <<<grid, block>>>(
                  Vx,Vz,hP,
                  xi_x, xi_z, zt_x, zt_z,
                  kappa3d, slw3d,
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
    sv_curv_col_ac_iso_rhs_cfspml(Vx,Vz,P,
                                  hVx,hVz,hP,
                                  xi_x, xi_z, zt_x, zt_z,
                                  kappa3d, slw3d,
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
      sv_curv_col_ac_iso_rhs_src_gpu  <<< grid,block >>> (
                  hVx,hVz,hP,
                  jac3d, slw3d, 
                  src_d,
                  verbose);
    }
  }
  // end func

  return 0;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

__global__ void
sv_curv_col_ac_iso_rhs_inner_gpu(
              float *  Vx , float *  Vz ,
              float *  P, 
              float * hVx , float * hVz ,
              float * hP, 
              float * xi_x, float * xi_z,
              float * zt_x, float * zt_z,
              float * kappa3d, float * slw3d,
              int ni1, int ni, int nk1, int nk,
              size_t siz_iz,
              int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
              int fdz_len, size_t * lfdz_shift, float * lfdz_coef,
              const int verbose)
{
  // local var
  float DxP,DxVx,DxVz;
  float DzP,DzVx,DzVz;
  float kappa,slw;
  float xix,xiz,ztx,ztz;

  float * Vx_ptr;
  float * Vz_ptr;
  float * P_ptr;

  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  size_t iz = blockIdx.y * blockDim.y + threadIdx.y;

  // caclu all points
  if(ix<ni && iz<nk)
  {
    size_t iptr = (ix+ni1) + (iz+nk1) * siz_iz;
    Vx_ptr = Vx + iptr;
    Vz_ptr = Vz + iptr;
    P_ptr  = P  + iptr;

    // Vx derivatives
    M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // Vz derivatives
    M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // P derivatives
    M_FD_SHIFT_PTR_MACDRP(DxP, P_ptr, fdx_len, lfdx_shift, lfdx_coef);
    M_FD_SHIFT_PTR_MACDRP(DzP, P_ptr, fdz_len, lfdz_shift, lfdz_coef);

    // metric
    xix = xi_x[iptr];
    xiz = xi_z[iptr];
    ztx = zt_x[iptr];
    ztz = zt_z[iptr];

    // medium
    kappa = kappa3d[iptr];
    slw = slw3d[iptr];

    // moment equation
    hVx[iptr] = - slw*( xix*DxP  
                       +ztx*DzP );
    hVz[iptr] = - slw*( xiz*DxP 
                       +ztz*DzP );

    // Hooke's equatoin
    hP[iptr] = -kappa *  ( xix*DxVx + ztx*DzVx
                          +xiz*DxVz + ztz*DzVz);
  }

  return;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

__global__ void
sv_curv_col_ac_iso_rhs_timg_z2_gpu(
               float *  P,
               int ni1, int ni, int nk1, int nk2, int nz,
               size_t siz_iz,
               const int verbose)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  for(int k=nk2; k<nz; k++)
  {
    if(ix<ni)
    {
      size_t iptr = (ix+ni1) + k * siz_iz;
      if(k == nk2)
      {
        P[iptr] = 0.0;
      }

      int k_phy = nk2 - (k-nk2);
      size_t iptr_phy = (ix+ni1) + k_phy * siz_iz;

      P[iptr] = -P[iptr_phy];
    }
  }

  return;
}

/*
 * implement vlow boundary
 */

__global__ void
sv_curv_col_ac_iso_rhs_vlow_z2_gpu(
               float *  Vx , float *  Vz ,
               float * hP, 
               float * xi_x, float * xi_z,
               float * zt_x, float * zt_z,
               float * kappa3d, float * slw3d,
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
  float kappa;
  float xix,xiz,ztx,ztz;

  float lfdz_coef[5] = {0.0};
  size_t   lfdz_shift[5] = {0};
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;

  // loop near surface layers
  //for (size_t n=0; n < 1; n++)
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
      kappa = kappa3d[iptr];

      // Vx derivatives
      M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      // Vz derivatives
      M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      if (k==nk2) // at surface, zero
      {
        hP[iptr] =  0.0;
      }
      else // lower than surface, lower order
      {
        M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
        M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);

        hP[iptr] = -kappa *  ( xix*DxVx + ztx*DzVx
                              +xiz*DxVz + ztz*DzVz);
      }
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
sv_curv_col_ac_iso_rhs_cfspml(
               float *  Vx , float *  Vz ,
               float *  P, 
               float * hVx , float * hVz ,
               float * hP,
               float * xi_x, float * xi_z,
               float * zt_x, float * zt_z,
               float * kappa3d, float * slw3d,
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

        sv_curv_col_ac_iso_rhs_cfspml_gpu <<<grid, block>>> (
                                idim, iside,
                                Vx, Vz, P, 
                                hVx, hVz, hP, 
                                xi_x, xi_z, zt_x, zt_z, 
                                kappa3d, slw3d,
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
sv_curv_col_ac_iso_rhs_cfspml_gpu(int idim, int iside,
                                  float * Vx,    float * Vz,
                                  float * P, 
                                  float * hVx,   float * hVz,
                                  float * hP, 
                                  float * xi_x,  float * xi_z,
                                  float * zt_x,  float * zt_z,
                                  float * kappa3d, float * slw3d,
                                  int nk2, size_t siz_iz,
                                  int fdx_len, size_t * lfdx_shift, float * lfdx_coef,
                                  int fdz_len, size_t * lfdz_shift, float * lfdz_coef,
                                  bdry_t bdry_d, const int verbose)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  size_t iz = blockIdx.y * blockDim.y + threadIdx.y;

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
  float DxP,DxVx,DxVz;
  float DzP,DzVx,DzVz;
  float kappa,slw;
  float xix,xiz,ztx,ztz;
  float hVx_rhs,hVz_rhs;
  float hP_rhs;

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
  float * pml_P    = abs_vars_cur + auxvar->Txx_pos;

  float * pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
  float * pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
  float * pml_hP   = abs_vars_rhs + auxvar->Txx_pos;

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
      kappa = kappa3d[iptr];
      slw = slw3d[iptr];

      // xi derivatives
      M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
      M_FD_SHIFT(DxP, P, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

      // combine for corr and aux vars
       hVx_rhs = -slw * ( xix*DxP                         );
       hVz_rhs = -slw * (                         xiz*DxP );
        hP_rhs = -kappa*( xix*DxVx            + xiz*DxVz  );

      // 1: make corr to moment equation
      hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
      hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

      // make corr to Hooke's equatoin
      hP[iptr] += coef_B_minus_1 * hP_rhs - coef_B * pml_P[iptr_a];
      
      // 2: aux var
      //   a1 = alpha + d / beta, dealt in abs_set_cfspml
      pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
      pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
      pml_hP[iptr_a] = coef_D * hP_rhs - coef_A * pml_P[iptr_a];
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
      kappa = kappa3d[iptr];

      // zt derivatives
      M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
      M_FD_SHIFT(DzP , P , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

      // combine for corr and aux vars
       hVx_rhs = -slw * ( ztx*DzP                         );
       hVz_rhs = -slw * (                         ztz*DzP );
        hP_rhs = -kappa*(ztx*DzVx            + ztz*DzVz);

      // 1: make corr to moment equation
      hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
      hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

      // make corr to Hooke's equatoin
      hP[iptr] += coef_B_minus_1 * hP_rhs - coef_B * pml_P[iptr_a];
      
      // 2: aux var
      //   a1 = alpha + d / beta, dealt in abs_set_cfspml
      pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
      pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
      pml_hP[iptr_a] = coef_D * hP_rhs - coef_A * pml_P[iptr_a];
    } 
  } // if which dim

  return;
}

/*******************************************************************************
 * add source terms
 ******************************************************************************/

__global__ void
sv_curv_col_ac_iso_rhs_src_gpu(
             float * hVx , float * hVz ,
             float * hP, 
             float * jac3d, float * slw3d,
             src_t src_d, // short nation for reference member
             const int verbose)
{
  // for easy coding and efficiency
  int max_ext = src_d.max_ext;

  // get fi / mij
  float fx, fz;
  float Mii;

  int it     = src_d.it;
  int istage = src_d.istage;
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;

  // add src; ix is a commont iterater var
  if(ix<src_d.total_number)
  {
    int   it_start = src_d.it_begin[ix];
    int   it_end   = src_d.it_end  [ix];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src_d.ext_indx + ix * max_ext;
      float *ptr_ext_coef = src_d.ext_coef + ix * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   ix * src_d.max_nt * src_d.max_stage // skip other src
                           + it_to_it_start * src_d.max_stage // skip other time step
                           + istage;
      if (src_d.force_actived == 1) {
        fx  = src_d.Fx [iptr_cur_stage];
        fz  = src_d.Fz [iptr_cur_stage];
      }
      if (src_d.moment_actived == 1) {
        Mii = src_d.Mxx[iptr_cur_stage];
      }
      
      // for extend points
      for (int i_ext=0; i_ext < src_d.ext_num[ix]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];

        if (src_d.force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          atomicAdd(&hVx[iptr], fx * V);
          atomicAdd(&hVz[iptr], fz * V);
        }

        if (src_d.moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          atomicAdd(&hP[iptr], -Mii * rjac);
        }
      } // i_ext

    } // it
  } 

  return;
}

