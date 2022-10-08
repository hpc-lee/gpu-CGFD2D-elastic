/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_ac_iso.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

int
sv_curv_col_ac_iso_onestage(
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
  float *__restrict__ P     = w_cur + wav->Txx_pos;
  float *__restrict__ hVx   = rhs   + wav->Vx_pos ; 
  float *__restrict__ hVz   = rhs   + wav->Vz_pos ; 
  float *__restrict__ hP    = rhs   + wav->Txx_pos; 

  float *__restrict__ xi_x  = metric->xi_x;
  float *__restrict__ xi_z  = metric->xi_z;
  float *__restrict__ zt_x  = metric->zeta_x;
  float *__restrict__ zt_z  = metric->zeta_z;
  float *__restrict__ jac3d = metric->jac;

  float *__restrict__ kappa3d = md->kappa;
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
  size_t siz_iz   = gd->siz_iz;

  // local fd op
  int              fdx_inn_len;
  int    *__restrict__ fdx_inn_indx;
  float  *__restrict__ fdx_inn_coef;
  int              fdz_inn_len;
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

  // free surface at z2 for pressure
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // imaging
    sv_curv_col_ac_iso_rhs_timg_z2(P,
                                   ni1,ni2,nk1,nk2,nz,
                                   siz_iz,
                                   verbose);
  }

  // inner points
  sv_curv_col_ac_iso_rhs_inner(Vx,Vz,P,
                               hVx,hVz,hP,
                               xi_x, xi_z, zt_x, zt_z,
                               kappa3d, slw3d,
                               ni1,ni2,nk1,nk2,siz_iz,
                               fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                               fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                               verbose);

  // free, abs, source in turn

  // free surface at z2 for velocity
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    // velocity: vlow
    sv_curv_col_ac_iso_rhs_vlow_z2(Vx,Vz,hP,
                                   xi_x, xi_z, zt_x, zt_z,
                                   kappa3d, slw3d,
                                   ni1,ni2,nk1,nk2,siz_iz,
                                   fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                   num_of_fdz_op,fdz_op,fdz_max_len,
                                   verbose);
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_curv_col_ac_iso_rhs_cfspml(Vx,Vz,P,
                                  hVx,hVz,hP,
                                  xi_x, xi_z, zt_x, zt_z,
                                  kappa3d, slw3d,
                                  nk2, siz_iz,
                                  fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                  fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                  bdry,
                                  verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_curv_col_ac_iso_rhs_src(hVx,hVz,hP,
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
sv_curv_col_ac_iso_rhs_inner(
              float *__restrict__  Vx , float *__restrict__  Vz ,
              float *__restrict__  P, 
              float *__restrict__ hVx , float *__restrict__ hVz ,
              float *__restrict__ hP, 
              float *__restrict__ xi_x, float *__restrict__ xi_z,
              float *__restrict__ zt_x, float *__restrict__ zt_z,
              float *__restrict__ kappa3d, float *__restrict__ slw3d,
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
  float DxP,DxVx,DxVz;
  float DzP,DzVx,DzVz;
  float kappa,slw;
  float xix,xiz,ztx,ztz;

  float *__restrict__ Vx_ptr;
  float *__restrict__ Vz_ptr;
  float *__restrict__ P_ptr;

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
        P_ptr  = P  + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // P derivatives
        M_FD_SHIFT_PTR_MACDRP(DxP, P_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzP, P_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

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

        iptr += 1;
      }
  }

  return 0;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

int
sv_curv_col_ac_iso_rhs_timg_z2(
               float *__restrict__  P,
               int ni1, int ni2, int nk1, int nk2, int nz,
               size_t siz_iz,
               const int verbose)
{
  // nk2
  size_t iptr_k = nk2 * siz_iz;
    size_t iptr = iptr_k + ni1;
    for (size_t i=ni1; i<=ni2; i++)
    {
      P[iptr] = 0.0;

      // next
      iptr += 1;
    }

  // mirror point
  for (size_t k=nk2+1; k<nz; k++)
  {
    int k_phy = nk2 - (k-nk2);
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr_gho = i + k     * siz_iz;
        size_t iptr_phy = i + k_phy * siz_iz;

        P[iptr_gho] = -P[iptr_phy];
      }
  }

  return 0;
}

/*
 * implement vlow boundary
 */

int
sv_curv_col_ac_iso_rhs_vlow_z2(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__ hP, 
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ kappa3d, float *__restrict__ slw3d,
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
  float kappa;
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
sv_curv_col_ac_iso_rhs_cfspml(
               float *__restrict__  Vx , float *__restrict__  Vz ,
               float *__restrict__  P, 
               float *__restrict__ hVx , float *__restrict__ hVz ,
               float *__restrict__ hP,
               float *__restrict__ xi_x, float *__restrict__ xi_z,
               float *__restrict__ zt_x, float *__restrict__ zt_z,
               float *__restrict__ kappa3d, float *__restrict__ slw3d,
               int nk2, size_t siz_iz,
               int fdx_len, int *__restrict__ fdx_indx, float *__restrict__ fdx_coef,
               int fdz_len, int *__restrict__ fdz_indx, float *__restrict__ fdz_coef,
               bdry_t *bdry,
               const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // val on point
  float DxP,DxVx,DxVz;
  float DzP,DzVx,DzVz;
  float kappa,slw;
  float xix,xiz,ztx,ztz;
  float hVx_rhs,hVz_rhs;
  float hP_rhs;

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
      float *__restrict__ pml_P    = abs_vars_cur + auxvar->Txx_pos;

      float *__restrict__ pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *__restrict__ pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *__restrict__ pml_hP   = abs_vars_rhs + auxvar->Txx_pos;

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
 * add source terms
 ******************************************************************************/

int
sv_curv_col_ac_iso_rhs_src(
             float *__restrict__ hVx , float *__restrict__ hVz ,
             float *__restrict__ hP, 
             float *__restrict__ jac3d, float *__restrict__ slw3d,
             src_t *src, // short nation for reference member
             const int verbose)
{
  int ierr = 0;

  // local var
  int si,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
  float fx, fz;
  float Mii;

  int it     = src->it;
  int istage = src->istage;

  // add src; is is a commont iterater var
  for (int is=0; is < src->total_number; is++)
  {
    int   it_start = src->it_begin[is];
    int   it_end   = src->it_end  [is];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src->ext_indx + is * max_ext;
      float *ptr_ext_coef = src->ext_coef + is * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   is * src->max_nt * src->max_stage // skip other src
                           + it_to_it_start * src->max_stage // skip other time step
                           + istage;
      if (src->force_actived == 1) {
        fx  = src->Fx [iptr_cur_stage];
        fz  = src->Fz [iptr_cur_stage];
      }
      if (src->moment_actived == 1) {
        Mii = src->Mxx[iptr_cur_stage];
      }
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          hVx[iptr] += fx * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          hP[iptr] -= Mii * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

