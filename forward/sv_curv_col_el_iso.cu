/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_el.h"
#include "sv_curv_col_el_iso.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

int
sv_curv_col_el_iso_onestage(
                  float *__restrict__ w_cur,
                  float *__restrict__ rhs, 
                  wav_t  *wav,
                  gd_t   *gd,
                  gdcurv_metric_t  *metric,
                  md_t *md,
                  bdry_t *bdry,
                  src_t *src,
                  // include different order/stentil
                  int num_of_fdx_op, fd_op_t *fdx_op,
                  int num_of_fdz_op, fd_op_t *fdz_op,
                  int fdz_max_len, 
                  const int verbose)
{
  // local pointer get each vars
  float *__restrict__ Vx    = w_cur + wav->Vx_pos ;
  float *__restrict__ Vz    = w_cur + wav->Vz_pos ;
  float *__restrict__ Txx   = w_cur + wav->Txx_pos;
  float *__restrict__ Tzz   = w_cur + wav->Tzz_pos;
  float *__restrict__ Txz   = w_cur + wav->Txz_pos;
  float *__restrict__ hVx   = rhs   + wav->Vx_pos ; 
  float *__restrict__ hVz   = rhs   + wav->Vz_pos ; 
  float *__restrict__ hTxx  = rhs   + wav->Txx_pos; 
  float *__restrict__ hTzz  = rhs   + wav->Tzz_pos; 
  float *__restrict__ hTxz  = rhs   + wav->Txz_pos; 

  float *__restrict__ xi_x  = metric->xi_x;
  float *__restrict__ xi_z  = metric->xi_z;
  float *__restrict__ zt_x  = metric->zeta_x;
  float *__restrict__ zt_z  = metric->zeta_z;
  float *__restrict__ jac3d = metric->jac;

  float *__restrict__ lam3d = md->lambda;
  float *__restrict__  mu3d = md->mu;
  float *__restrict__ slw3d = md->rho;

  // grid size
  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  int nk1 = gd->nk1;
  int nk2 = gd->nk2;

  int ni  = gd->ni;
  int nk  = gd->nk;
  int nx  = gd->nx;
  int nz  = gd->nz;
  size_t siz_iz  = gd->siz_iz;

  float *vecVx2Vz = bdry->vecVx2Vz2;

  // local fd op
  int                  fdx_inn_len;
  int    *__restrict__ fdx_inn_indx;
  float  *__restrict__ fdx_inn_coef;
  int                  fdz_inn_len;
  int    *__restrict__ fdz_inn_indx;
  float  *__restrict__ fdz_inn_coef;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op
  fdx_inn_len  = fdx_op[num_of_fdx_op-1].total_len;
  fdx_inn_indx = fdx_op[num_of_fdx_op-1].indx;
  fdx_inn_coef = fdx_op[num_of_fdx_op-1].coef;

  fdz_inn_len  = fdz_op[num_of_fdz_op-1].total_len;
  fdz_inn_indx = fdz_op[num_of_fdz_op-1].indx;
  fdz_inn_coef = fdz_op[num_of_fdz_op-1].coef;

  // inner points
  sv_curv_col_el_iso_rhs_inner(Vx,Vz,Txx,Tzz,Txz,
                               hVx,hVz,hTxx,hTzz,hTxz,
                               xi_x, xi_z, zt_x, zt_z,
                               lam3d, mu3d, slw3d,
                               ni1,ni2,nk1,nk2,siz_iz,
                               fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                               fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                               verbose);

  // free, abs, source in turn

  // free surface at z2
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // tractiong
    sv_curv_col_el_rhs_timg_z2(Txx,Tzz,Txz,hVx,hVz,
                               xi_x, xi_z, zt_x, zt_z,
                               jac3d, slw3d,
                               ni1,ni2,nk1,nk2,siz_iz,
                               fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                               fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                               verbose);

    // velocity: vlow
    sv_curv_col_el_iso_rhs_vlow_z2(Vx,Vz,hTxx,hTzz,hTxz,
                                   xi_x, xi_z, zt_x, zt_z,
                                   lam3d, mu3d, slw3d,
                                   vecVx2Vz,
                                   ni1,ni2,nk1,nk2,siz_iz,
                                   fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                   num_of_fdz_op,fdz_op,fdz_max_len,
                                   verbose);
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_curv_col_el_iso_rhs_cfspml(Vx,Vz,Txx,Tzz,Txz,
                                  hVx,hVz,hTxx,hTzz,hTxz,
                                  xi_x, xi_z, zt_x, zt_z,
                                  lam3d, mu3d, slw3d,
                                  nk2, siz_iz,
                                  fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                  fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                  bdry,
                                  verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_curv_col_el_rhs_src(hVx,hVz,hTxx,hTzz,hTxz,
                           jac3d, slw3d, 
                           src,
                           verbose);
  }
  // end func

  return 0;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

int
sv_curv_col_el_iso_rhs_inner(
                float *__restrict__  Vx , float *__restrict__  Vz ,
                float *__restrict__  Txx, float *__restrict__  Tzz,
                float *__restrict__  Txz, 
                float *__restrict__ hVx , float *__restrict__ hVz ,
                float *__restrict__ hTxx, float *__restrict__ hTzz,
                float *__restrict__ hTxz, 
                float *__restrict__ xi_x, float *__restrict__ xi_z,
                float *__restrict__ zt_x, float *__restrict__ zt_z,
                float *__restrict__ lam3d, float *__restrict__ mu3d, float *__restrict__ slw3d,
                int ni1, int ni2, int nk1, int nk2,
                size_t siz_iz, 
                int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
                int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
                const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;

  float *__restrict__ Vx_ptr;
  float *__restrict__ Vz_ptr;
  float *__restrict__ Txx_ptr;
  float *__restrict__ Txz_ptr;
  float *__restrict__ Tzz_ptr;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
        Vx_ptr = Vx + iptr;
        Vz_ptr = Vz + iptr;
        Txx_ptr = Txx + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxx, Txx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTzz, Tzz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // metric
        xix = xi_x[iptr];
        xiz = xi_z[iptr];
        ztx = zt_x[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // moment equation
        hVx[iptr] = slw*( xix*DxTxx + xiz*DxTxz  
                         +ztx*DzTxx + ztz*DzTxz );
        hVz[iptr] = slw*( xix*DxTxz + xiz*DxTzz 
                         +ztx*DzTxz + ztz*DzTzz );

        // Hooke's equatoin
        hTxx[iptr] =  lam2mu * ( xix*DxVx + ztx*DzVx)
                    + lam    * ( xiz*DxVz + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz  + ztz*DzVz)
                    +lam    * ( xix*DxVx  + ztx*DzVx);

        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +ztz*DzVx + ztx*DzVz
                    );

        iptr += 1;
      }
  }

  return 0;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement vlow boundary
 */

int
sv_curv_col_el_iso_rhs_vlow_z2(
                   float *__restrict__  Vx , float *__restrict__  Vz ,
                   float *__restrict__ hTxx, float *__restrict__ hTzz,
                   float *__restrict__ hTxz, 
                   float *__restrict__ xi_x, float *__restrict__ xi_z,
                   float *__restrict__ zt_x, float *__restrict__ zt_z,
                   float *__restrict__ lam3d, float *__restrict__ mu3d, float *__restrict__ slw3d,
                   float *__restrict__ vecVx2Vz, 
                   int ni1, int ni2, int nk1, int nk2,
                   size_t siz_iz, 
                   int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
                   int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
                   const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  int    lfdz_shift[fdz_max_len];

  // local var
  int i,k;
  int n_fd; // loop var for fd
  int fdz_len;

  // local var
  float DxVx,DxVz;
  float DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }

  // loop near surface layers
  //for (size_t n=0; n < 1; n++)
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    int  lfdz_len  = fdz_op[n].total_len;
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_op[n].indx;
    float *p_fdz_coef  = fdz_op[n].coef;
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = p_fdz_indx[n_fd] * siz_iz;
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd];
    }

    // for index
    size_t iptr_k = k * siz_iz;

      size_t iptr = iptr_k + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // metric
        xix = xi_x[iptr];
        xiz = xi_z[iptr];
        ztx = zt_x[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        if (k==nk2) // at surface, convert
        {
          size_t ij = (i )*4;
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
        hTxx[iptr] =  lam2mu * ( xix*DxVx  + ztx*DzVx)
                    + lam    * ( xiz*DxVz  + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz  + ztz*DzVz)
                    +lam    * ( xix*DxVx  + ztx*DzVx);

        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +ztz*DzVx + ztx*DzVz
                    );

        iptr += 1;
      }
  }

  return 0;
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

int
sv_curv_col_el_iso_rhs_cfspml(
                   float *__restrict__  Vx , float *__restrict__  Vz ,
                   float *__restrict__  Txx, float *__restrict__  Tzz,
                   float *__restrict__  Txz, 
                   float *__restrict__ hVx , float *__restrict__ hVz ,
                   float *__restrict__ hTxx, float *__restrict__ hTzz,
                   float *__restrict__ hTxz, 
                   float *__restrict__ xi_x, float *__restrict__ xi_z,
                   float *__restrict__ zt_x, float *__restrict__ zt_z,
                   float *__restrict__ lam3d, float *__restrict__  mu3d, float *__restrict__ slw3d,
                   int nk2, size_t siz_iz, 
                   int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
                   int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
                   bdry_t *bdry,
                   const int verbose)
{
  float *vecVx2Vz = bdry->vecVx2Vz2;

  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTzz,DxTxz,DxVx,DxVz;
  float DzTxx,DzTzz,DzTxz,DzVx,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiz,ztx,ztz;
  float hVx_rhs,hVz_rhs;
  float hTxx_rhs,hTzz_rhs,hTxz_rhs;
  // for free surface
  float Dx_DzVx,Dx_DzVz;

  // local
  int i,k;
  int iptr, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_B_minus_1;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_iz;
  }

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip to next face if not cfspml
      if (bdry->is_sides_pml[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdry->ni1[idim][iside];
      int abs_ni2 = bdry->ni2[idim][iside];
      int abs_nk1 = bdry->nk1[idim][iside];
      int abs_nk2 = bdry->nk2[idim][iside];

      // get coef for this face
      float *__restrict__ ptr_coef_A = bdry->A[idim][iside];
      float *__restrict__ ptr_coef_B = bdry->B[idim][iside];
      float *__restrict__ ptr_coef_D = bdry->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);

      // get pml vars
      float *__restrict__ abs_vars_cur = auxvar->cur;
      float *__restrict__ abs_vars_rhs = auxvar->rhs;

      float *__restrict__ pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
      float *__restrict__ pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
      float *__restrict__ pml_Txx  = abs_vars_cur + auxvar->Txx_pos;
      float *__restrict__ pml_Tzz  = abs_vars_cur + auxvar->Tzz_pos;
      float *__restrict__ pml_Txz  = abs_vars_cur + auxvar->Txz_pos;

      float *__restrict__ pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *__restrict__ pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *__restrict__ pml_hTxx = abs_vars_rhs + auxvar->Txx_pos;
      float *__restrict__ pml_hTzz = abs_vars_rhs + auxvar->Tzz_pos;
      float *__restrict__ pml_hTxz = abs_vars_rhs + auxvar->Txz_pos;

      // for each dim
      if (idim == 0 ) // x direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;
            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // pml coefs
              int abs_i = i - abs_ni1;
              coef_D = ptr_coef_D[abs_i];
              coef_A = ptr_coef_A[abs_i];
              coef_B = ptr_coef_B[abs_i];
              coef_B_minus_1 = coef_B - 1.0;

              // metric
              xix = xi_x[iptr];
              xiz = xi_z[iptr];
              //xiz = fabs(xiz);

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              slw = slw3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // xi derivatives
              M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( xix*DxTxx + xiz*DxTxz );
               hVz_rhs = slw * ( xix*DxTxz + xiz*DxTzz );
              hTxx_rhs = lam2mu*xix*DxVx + lam*xiz*DxVz;
              hTzz_rhs = lam*xix*DxVx + lam2mu*xiz*DxVz;
              hTxz_rhs = mu*( xiz*DxVx + xix*DxVz );

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
              if (bdry->is_sides_free[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i )*4;
                Dx_DzVx = vecVx2Vz[ij+2*0+0] * DxVx
                        + vecVx2Vz[ij+2*0+1] * DxVz;

                Dx_DzVz = vecVx2Vz[ij+2*1+0] * DxVx
                        + vecVx2Vz[ij+2*1+1] * DxVz;

                // metric
                ztx = zt_x[iptr];
                ztz = zt_z[iptr];

                // keep xi derivative terms, including free surface convered
                hTxx_rhs =    lam2mu * (            ztx*Dx_DzVx)
                            + lam    * (            ztz*Dx_DzVz);

                hTzz_rhs =   lam2mu * (            ztz*Dx_DzVz)
                            +lam    * (            ztx*Dx_DzVx);

                hTxz_rhs = mu *(
                             ztz*Dx_DzVx + ztx*Dx_DzVz
                            );

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
              }

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      }
      else // z direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_iz;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_B_minus_1 = coef_B - 1.0;

            iptr = iptr_k + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // metric
              ztx = zt_x[iptr];
              ztz = zt_z[iptr];
              //ztx = fabs(ztx);

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              slw = slw3d[iptr];
              lam2mu = lam + 2.0 * mu;

              // zt derivatives
              M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( ztx*DzTxx + ztz*DzTxz );
               hVz_rhs = slw * ( ztx*DzTxz + ztz*DzTzz );
              hTxx_rhs = lam2mu*ztx*DzVx + lam*ztz*DzVz;
              hTzz_rhs = lam*ztx*DzVx + lam2mu*ztz*DzVz;
              hTxz_rhs = mu*( ztz*DzVx + ztx*DzVz );

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

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
        } // k
      } // if which dim
    } // iside
  } // idim

  return 0;
}

/*******************************************************************************
 * free surface coef
 ******************************************************************************/

int
sv_curv_col_el_iso_dvh2dvz(gd_t        *gd,
                           gdcurv_metric_t *metric,
                           md_t       *md,
                           bdry_t     *bdry,
                           const int verbose)
{
  int ierr = 0;

  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  int nk1 = gd->nk1;
  int nk2 = gd->nk2;
  int nx  = gd->nx;
  int nz  = gd->nz;
  size_t siz_iz  = gd->siz_iz;

  // point to each var
  float *__restrict__ xi_x = metric->xi_x;
  float *__restrict__ xi_z = metric->xi_z;
  float *__restrict__ zt_x = metric->zeta_x;
  float *__restrict__ zt_z = metric->zeta_z;

  float *__restrict__ lam3d = md->lambda;
  float *__restrict__  mu3d = md->mu;

  float *vecVx2Vz = bdry->vecVx2Vz2;
  
  float A[2][2], B[2][2];
  float AB[2][2];

  int k = nk2;

    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + k * siz_iz;

      float e11 = xi_x[iptr];
      float e12 = xi_z[iptr];
      float e21 = zt_x[iptr];
      float e22 = zt_z[iptr];

      float lam    = lam3d[iptr];
      float miu    =  mu3d[iptr];
      float lam2mu = lam + 2.0f * miu;

      // first dim: irow; sec dim: jcol, as Fortran code
      A[0][0]=lam2mu*e21*e21+miu*e22*e22;
      A[0][1]=lam*e21*e22+miu*e22*e21;
      A[1][0]=lam*e22*e21+miu*e21*e22;
      A[1][1]=lam2mu*e22*e22+miu*e21*e21;

      fdlib_math_invert2x2(A);

      B[0][0]=-lam2mu*e21*e11-miu*e22*e12;
      B[0][1]=-lam*e21*e12-miu*e22*e11;
      B[1][0]=-lam*e22*e11-miu*e21*e12;
      B[1][1]=-lam2mu*e22*e12-miu*e21*e11;

      fdlib_math_matmul2x2(A, B, AB);

      size_t ij = (i) * 4;

      // save into mat
      for(int irow = 0; irow < 2; irow++)
        for(int jcol = 0; jcol < 2; jcol++){
          vecVx2Vz[ij + irow*2 + jcol] = AB[irow][jcol];
        }
    }

  return ierr;
}

