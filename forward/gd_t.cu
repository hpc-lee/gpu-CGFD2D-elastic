/*
********************************************************************************
* Curve grid metric calculation using MacCormack scheme                        *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "fd_t.h"
#include "gd_t.h"

// used in read grid file
#define M_gd_INDEX( i, k, ni ) ( ( i ) + ( k ) * ( ni ) )

void 
gd_curv_init(gd_t *gdcurv)
{
  /*
   * 0-2: x2d, z2d
   */

  gdcurv->type = GD_TYPE_CURV;

  gdcurv->ncmp = CONST_NDIM;

  gdcurv->siz_iz   = gdcurv->nx;
  gdcurv->siz_icmp = gdcurv->nx * gdcurv->nz;
  
  // vars
  gdcurv->v3d = (float *) fdlib_mem_calloc_1d_float(
                  gdcurv->siz_icmp * gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // position of each v3d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // name of each v3d
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(gdcurv->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_init");
  
  // set value
  int icmp = 0;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","x");
  gdcurv->x2d = gdcurv->v3d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","z");
  gdcurv->z2d = gdcurv->v3d + cmp_pos[icmp];
  
  // set pointer
  gdcurv->cmp_pos  = cmp_pos;
  gdcurv->cmp_name = cmp_name;

  return;
}

void 
gd_curv_metric_init(gd_t        *gdcurv,
                    gdcurv_metric_t *metric)
{
  const int num_grid_vars = 5; 
  /*
   * 0: jac
   * 1-2: xi_x,  xi_z
   * 3-4: zeta_x, zeta_z
   */

  metric->nx   = gdcurv->nx;
  metric->nz   = gdcurv->nz;
  metric->ncmp = num_grid_vars;

  metric->siz_iz   = metric->nx;
  metric->siz_icmp  = metric->nx * metric->nz;
  
  // vars
  metric->v3d = (float *) fdlib_mem_calloc_1d_float(
                  metric->siz_icmp * metric->ncmp, 0.0, "gd_curv_init_g3d");
  if (metric->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc metric vars\n");
      fflush(stderr);
  }
  
  // position of each v3d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(metric->ncmp,
                                                         0, 
                                                         "gd_curv_metric_init");
  
  // name of each v3d
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(metric->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_metric_init");
  
  // set value
  for (int icmp=0; icmp < metric->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * metric->siz_icmp;
  }

  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","jac");
  metric->jac = metric->v3d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","xi_x");
  metric->xi_x = metric->v3d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","xi_z");
  metric->xi_z = metric->v3d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","zeta_x");
  metric->zeta_x = metric->v3d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","zeta_z");
  metric->zeta_z = metric->v3d + cmp_pos[icmp];
  
  // set pointer
  metric->cmp_pos  = cmp_pos;
  metric->cmp_name = cmp_name;
}

//
// need to change to use fdlib_math.c
//
void
gd_curv_metric_cal(gd_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int * fd_indx, float * fd_coef)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx  = gdcurv->nx;
  int nz  = gdcurv->nz;
  size_t siz_iz  = gdcurv->siz_iz;

  // point to each var
  float * x2d  = gdcurv->x2d;
  float * z2d  = gdcurv->z2d;
  float * jac2d= metric->jac;
  float * xi_x = metric->xi_x;
  float * xi_z = metric->xi_z;
  float * zt_x = metric->zeta_x;
  float * zt_z = metric->zeta_z;

  float x_xi, x_zt;
  float z_xi, z_zt;
  float jac;
  float vec1[3], vec2[3], vec3[3], vecg[3];
  int n_fd;

  // use local stack array for speedup
  float  lfd_coef [fd_len];
  int    lfdx_shift[fd_len];
  int    lfdz_shift[fd_len];
  // put fd op into local array
  for (int k=0; k < fd_len; k++) {
    lfd_coef [k] = fd_coef[k];
    lfdx_shift[k] = fd_indx[k]            ;
    lfdz_shift[k] = fd_indx[k] * siz_iz;
  }

  for (size_t k = nk1; k <= nk2; k++){
      for (size_t i = ni1; i <= ni2; i++)
      {
        size_t iptr = i + k * siz_iz;

        x_xi = 0.0; x_zt = 0.0;
        z_xi = 0.0; z_zt = 0.0;

        M_FD_SHIFT(x_xi, x2d, iptr, fd_len, lfdx_shift, lfd_coef, n_fd);
        M_FD_SHIFT(z_xi, z2d, iptr, fd_len, lfdx_shift, lfd_coef, n_fd);

        M_FD_SHIFT(x_zt, x2d, iptr, fd_len, lfdz_shift, lfd_coef, n_fd);
        M_FD_SHIFT(z_zt, z2d, iptr, fd_len, lfdz_shift, lfd_coef, n_fd);

        vec1[0] = x_xi; vec1[1] = z_xi; vec1[2] = 0.0;
        vec2[0] = x_zt; vec2[1] = z_zt; vec2[2] = 0.0;
        vec3[0] = 0.0 ; vec3[1] = 0.0 ; vec3[2] = 1.0;

        // jac
        fdlib_math_cross_product(vec1, vec2, vecg);
        jac = fdlib_math_dot_product(vecg, vecg);
        jac = sqrt(jac);
        jac2d[iptr]  = jac;

        // xi_i
        fdlib_math_cross_product(vec2, vec3, vecg);
        xi_x[iptr] = vecg[0] / jac;
        xi_z[iptr] = vecg[1] / jac;

        // zt_i
        fdlib_math_cross_product(vec3, vec1, vecg);
        zt_x[iptr] = vecg[0] / jac;
        zt_z[iptr] = vecg[1] / jac;
      }
  }

  mirror_symmetry(gdcurv,metric->v3d,metric->ncmp);
  //geometric_symmetry(gdcurv,metric->v3d,metric->ncmp);

  return;
}

int mirror_symmetry(gd_t *gdcurv,float *v3d, int ncmp)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx  = gdcurv->nx;
  int nz  = gdcurv->nz;
  size_t siz_iz  = gdcurv->siz_iz;
  size_t siz_icmp  = gdcurv->siz_icmp;

  size_t iptr, iptr1, iptr2; 
  for(int icmp=0; icmp<ncmp; icmp++){
    iptr = icmp * siz_icmp;
    // x1, mirror
    for (size_t k = 0; k < nz; k++){
      for (size_t i = 0; i < ni1; i++)
      {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + (2*ni1-i) + k * siz_iz;
        v3d[iptr1] = v3d[iptr2];
      }
    }
    // x2, mirror
    for (size_t k = 0; k < nz; k++){
      for (size_t i = ni2+1; i < nx; i++)
      {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + (2*ni2-i) + k * siz_iz;
        v3d[iptr1] = v3d[iptr2];
      }
    }
    // z1, mirror
    for (size_t k = 0; k < nk1; k++) {
      for (size_t i = 0; i < nx; i++) {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + i + (2*nk1-k) * siz_iz;
        v3d[iptr1] = v3d[iptr2];
      }
    }
    // z2, mirror
    for (size_t k = nk2+1; k < nz; k++) {
      for (size_t i = 0; i < nx; i++) {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + i + (2*nk2-k) * siz_iz;
        v3d[iptr1] = v3d[iptr2];
      }
    }
  }

  return 0;
}

int geometric_symmetry(gd_t *gdcurv,float *v3d, int ncmp)
{
  int ni1 = gdcurv->ni1;
  int ni2 = gdcurv->ni2;
  int nk1 = gdcurv->nk1;
  int nk2 = gdcurv->nk2;
  int nx  = gdcurv->nx;
  int nz  = gdcurv->nz;
  size_t siz_iz  = gdcurv->siz_iz;
  size_t siz_icmp  = gdcurv->siz_icmp;

  size_t iptr, iptr1, iptr2, iptr3; 
  for(int icmp=0; icmp<ncmp; icmp++){
    iptr = icmp * siz_icmp;
    // x1 
    for (size_t k = 0; k < nz; k++){
      for (size_t i = 0; i < ni1; i++)
      {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + ni1 + k * siz_iz;
        iptr3 = iptr + (2*ni1-i) + k * siz_iz;
        v3d[iptr1] = 2*v3d[iptr2] - v3d[iptr3];
      }
    }
    // x2, mirror
    for (size_t k = 0; k < nz; k++){
      for (size_t i = ni2+1; i < nx; i++)
      {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + ni2 + k * siz_iz;
        iptr3 = iptr + (2*ni2-i) + k * siz_iz;
        v3d[iptr1] = 2*v3d[iptr2] - v3d[iptr3];
      }
    }
    // z1, mirror
    for (size_t k = 0; k < nk1; k++) {
      for (size_t i = 0; i < nx; i++) {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + i + nk1 * siz_iz;
        iptr3 = iptr + i + (2*nk1-k) * siz_iz;
        v3d[iptr1] = 2*v3d[iptr2] - v3d[iptr3];
      }
    }
    // z2, mirror
    for (size_t k = nk2+1; k < nz; k++) {
      for (size_t i = 0; i < nx; i++) {
        iptr1 = iptr + i + k * siz_iz;
        iptr2 = iptr + i + nk2 * siz_iz;
        iptr3 = iptr + i + (2*nk2-k) * siz_iz;
        v3d[iptr1] = 2*v3d[iptr2] - v3d[iptr3];
      }
    }
  }

  return 0;
}

/*
 * generate cartesian grid for curv struct
 */
void
gd_curv_gen_cart(gd_t *gdcurv,
                 float dx, float x0_glob,
                 float dz, float z0_glob)
{
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;

  float x0 = x0_glob + (0 - gdcurv->fdx_nghosts) * dx;
  float z0 = z0_glob + (0 - gdcurv->fdz_nghosts) * dz;

  size_t iptr = 0;
  for (size_t k=0; k<gdcurv->nz; k++)
  {
      for (size_t i=0; i<gdcurv->nx; i++)
      {
        x2d[iptr] = x0 + i * dx;
        z2d[iptr] = z0 + k * dz;

        iptr++;
      }
  }

  return;
}

/*
 * generate cartesian grid
 */

void 
gd_cart_init_set(gd_t *gdcart,
                 float dx, float x0_glob,
                 float dz, float z0_glob)
{
  /*
   * 0-2: x2d, z2d
   */

  gdcart->type = GD_TYPE_CART;

  gdcart->nx   = gdcart->nx;
  gdcart->nz   = gdcart->nz;
  gdcart->ncmp = CONST_NDIM;

  gdcart->siz_iz   = gdcart->nx;
  gdcart->siz_icmp = gdcart->nx * gdcart->nz;
  
  // vars
  float *x1d = (float *) fdlib_mem_calloc_1d_float(
                  gdcart->nx, 0.0, "gd_cart_init");
  float *z1d = (float *) fdlib_mem_calloc_1d_float(
                  gdcart->nz, 0.0, "gd_cart_init");
  if (z1d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  float x0 = x0_glob + (0 - gdcart->fdx_nghosts) * dx;
  float z0 = z0_glob + (0 - gdcart->fdz_nghosts) * dz;

  for (size_t k=0; k< gdcart->nz; k++)
  {
        z1d[k] = z0 + k * dz;
  }
  for (size_t i=0; i< gdcart->nx; i++)
  {
        x1d[i] = x0 + i * dx;
  }

  gdcart->dx = dx;
  gdcart->dz = dz;

  gdcart->xmin = x0;
  gdcart->zmin = z0;
  gdcart->xmax = x0 + (gdcart->nx-1) * dx;
  gdcart->zmax = z0 + (gdcart->nz-1) * dz;

  gdcart->x0_glob = x0_glob;
  gdcart->z0_glob = z0_glob;

  gdcart->x1d = x1d;
  gdcart->z1d = z1d;
  
  return;
}

//
// input/output
//
void
gd_curv_coord_export(gd_t *gdcurv,
                     char *output_dir)
{
  int ni  = gdcurv->ni;
  int nk  = gdcurv->nk;
  int ni1 = gdcurv->ni1;
  int nk1 = gdcurv->nk1;
  int ni2 = gdcurv->ni2;
  int nk2 = gdcurv->nk2;
  size_t iptr, iptr1;
  size_t siz_iz = gdcurv->siz_iz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  float *coord_x = (float *) malloc(sizeof(float)*ni*nk);  
  float *coord_z = (float *) malloc(sizeof(float)*ni*nk);  

  for (int k=nk1; k<=nk2; k++) {
    for (int i=ni1; i<=ni2; i++) {
      iptr = i + siz_iz * k;
      iptr1 = (i-3) + ni * (k-3);
      coord_x[iptr1] = x2d[iptr]; 
      coord_z[iptr1] = z2d[iptr];
    }
  }
  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  int ncid;
  int xid,zid;
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", ni, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nk, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, CONST_NDIM, dimid, &xid);
  handle_nc_err(ierr);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, CONST_NDIM, dimid, &zid);
  handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, xid, coord_x);
  handle_nc_err(ierr);
  ierr = nc_put_var_float(ncid, zid, coord_z);
  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  free(coord_x);
  free(coord_z);

  return;
}

void
gd_curv_coord_import(gd_t *gdcurv, char *import_dir)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/coord.nc", import_dir);
  
  int ni = gdcurv->ni;
  int nk = gdcurv->nk;
  int ni1 = gdcurv->ni1;
  int nk1 = gdcurv->nk1;
  int ni2 = gdcurv->ni2;
  int nk2 = gdcurv->nk2;
  size_t siz_iz = gdcurv->siz_iz;
  float *x2d = gdcurv->x2d;
  float *z2d = gdcurv->z2d;
  
  float *coord_x = (float *) malloc(sizeof(float)*ni*nk);  
  float *coord_z = (float *) malloc(sizeof(float)*ni*nk);  
  size_t start[] = {0, 0};
  size_t count[] = {nk, ni};
    
  // read grid nc file
  int ncid;
  int xid,zid;
  int ierr, iptr, iptr1;
   
  ierr = nc_open(in_file, NC_NOWRITE, &ncid);  handle_nc_err(ierr);
  
  // read vars
  ierr = nc_inq_varid(ncid, "x", &xid);  handle_nc_err(ierr);
  ierr = nc_inq_varid(ncid, "z", &zid);  handle_nc_err(ierr);
  
  ierr = nc_get_vara_float(ncid, xid, start, count, coord_x); handle_nc_err(ierr);
  ierr = nc_get_vara_float(ncid, zid, start, count, coord_z); handle_nc_err(ierr);
                                                           
  // close file
  ierr = nc_close(ncid); handle_nc_err(ierr);
  
  for (int k = nk1; k <= nk2; k++){
    for (int i = ni1; i<= ni2; i++){
      iptr = i + siz_iz * k;
      iptr1 = (i-3) + ni * (k-3);
      x2d[iptr] = coord_x[iptr1]; 
      z2d[iptr] = coord_z[iptr1];
    }
  }
  
  geometric_symmetry(gdcurv,gdcurv->v3d,gdcurv->ncmp);
  
  free(coord_x);
  free(coord_z);

  return;
}

/*
void
gd_curv_coord_export(gd_t *gdcurv,
                     char *output_dir)
{
  size_t * c3d_pos   = gdcurv->cmp_pos;
  char  ** c3d_name  = gdcurv->cmp_name;
  int number_of_vars = gdcurv->ncmp;
  int  nx = gdcurv->nx;
  int  nz = gdcurv->nz;
  int  ni1 = gdcurv->ni1;
  int  nk1 = gdcurv->nk1;
  int  ni  = gdcurv->ni;
  int  nk  = gdcurv->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[gdcurv->ncmp];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v3d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return;
}

void
gd_curv_coord_import(gd_t *gdcurv, char *import_dir)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/coord.nc", import_dir);
  
  // read in nc
  int ncid;
  int varid;

  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);  handle_nc_err(ierr);

  // read vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++)
  {
    float *ptr = gdcurv->v3d + gdcurv->cmp_pos[ivar];

    ierr = nc_inq_varid(ncid, gdcurv->cmp_name[ivar], &varid);  handle_nc_err(ierr);

    ierr = nc_get_var(ncid, varid, ptr);  handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);  handle_nc_err(ierr);

  return;
}
*/

void
gd_cart_coord_export(gd_t *gdcart,
                     char *output_dir)
{
  int  nx = gdcart->nx;
  int  nz = gdcart->nz;
  int  ni1 = gdcart->ni1;
  int  nk1 = gdcart->nk1;
  int  ni  = gdcart->ni;
  int  nk  = gdcart->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[CONST_NDIM];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, 1, dimid+1, &varid[0]);
  handle_nc_err(ierr);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, 1, dimid+0, &varid[1]);
  handle_nc_err(ierr);

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  ierr = nc_put_var_float(ncid, varid[0], gdcart->x1d);
  handle_nc_err(ierr);
  ierr = nc_put_var_float(ncid, varid[1], gdcart->z1d);
  handle_nc_err(ierr);
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return;
}

void
gd_curv_metric_export(gd_t        *gd,
                      gdcurv_metric_t *metric,
                      char *output_dir)
{
  size_t * g3d_pos   = metric->cmp_pos;
  char  ** g3d_name  = metric->cmp_name;
  int  number_of_vars = metric->ncmp;
  int  nx = metric->nx;
  int  nz = metric->nz;
  int  ni1 = gd->ni1;
  int  nk1 = gd->nk1;
  int  ni  = gd->ni;
  int  nk  = gd->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/metric.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, g3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = metric->v3d + g3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);
}

void
gd_curv_metric_import(gdcurv_metric_t *metric, char *import_dir)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/metric.nc", import_dir);
  
  // read in nc
  int ncid;
  int varid;

  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  handle_nc_err(ierr);

  // read vars
  for (int ivar=0; ivar<metric->ncmp; ivar++)
  {
    float *ptr = metric->v3d + metric->cmp_pos[ivar];

    ierr = nc_inq_varid(ncid, metric->cmp_name[ivar], &varid);
    handle_nc_err(ierr);

    ierr = nc_get_var(ncid, varid, ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return;
}

void
gd_curv_set_minmax(gd_t *gdcurv)
{
  float xmin = gdcurv->x2d[0], xmax = gdcurv->x2d[0];
  float zmin = gdcurv->z2d[0], zmax = gdcurv->z2d[0];
  
  for (size_t i = 0; i < gdcurv->siz_icmp; i++){
      xmin = xmin < gdcurv->x2d[i] ? xmin : gdcurv->x2d[i];
      xmax = xmax > gdcurv->x2d[i] ? xmax : gdcurv->x2d[i];
      zmin = zmin < gdcurv->z2d[i] ? zmin : gdcurv->z2d[i];
      zmax = zmax > gdcurv->z2d[i] ? zmax : gdcurv->z2d[i];
  }

  gdcurv->xmin = xmin;
  gdcurv->xmax = xmax;
  gdcurv->zmin = zmin;
  gdcurv->zmax = zmax;

  return;
}

/*
 * convert cart coord to global index
 */

int
gd_cart_coord_to_local_indx(gd_t *gdcart,
                            float sx,
                            float sz,
                            int   *ou_si, int *ou_sk,
                            float *ou_sx_inc, float *ou_sz_inc)
{
  int ierr = 0;

  int si_glob = (int)( (sx - gdcart->x0_glob) / gdcart->dx + 0.5 );
  int sk_glob = (int)( (sz - gdcart->z0_glob) / gdcart->dz + 0.5 );
  float sx_inc = si_glob * gdcart->dx + gdcart->x0_glob - sx;
  float sz_inc = sk_glob * gdcart->dz + gdcart->z0_glob - sz;

  *ou_si = si_glob + gdcart->fdx_nghosts;
  *ou_sk = sk_glob + gdcart->fdz_nghosts;
  *ou_sx_inc = sx_inc;
  *ou_sz_inc = sz_inc;

  return ierr; 
}

/* 
 * if the nearest point in this thread then search its grid index
 *   return value:
 *      1 - in this thread
 *      0 - not in this thread
 */

int
gd_curv_coord_to_local_indx(gd_t *gd,
                            float sx, float sz,
                            int *si, int *sk,
                            float *sx_inc, float *sz_inc) 
{
  int is_here = 0; // default outside

  int nx = gd->nx;
  int nz = gd->nz;
  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  int nk1 = gd->nk1;
  int nk2 = gd->nk2;
  size_t siz_iz = gd->siz_iz;

  float * x2d = gd->x2d;
  float * z2d = gd->z2d;

  // outside coord range
  if ( sx < gd->xmin || sx > gd->xmax ||
       sz < gd->zmin || sz > gd->zmax)
  {
    is_here = 0;
    *si = -1000;
    *sk = -1000;
    *sx_inc = 0.0;
    *sz_inc = 0.0;
    return is_here;
  }

  // init closest point
  float min_dist = sqrtf(  (sx - x2d[0]) * (sx - x2d[0])
      + (sz - z2d[0]) * (sz - z2d[0]) );
  int min_dist_i = 0 ;
  int min_dist_k = 0 ;

  // compute distance to each point
  for (int k=0; k<nz; k++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float x = x2d[iptr];
        float z = z2d[iptr];

        float DistInt = sqrtf(  (sx - x) * (sx - x)
            + (sz - z) * (sz - z) );

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = i;
          min_dist_k = k;
        }
      }
  }

  // if nearest index is outside phys region, not here
  if ( min_dist_i < ni1 || min_dist_i > ni2 ||
      min_dist_k < nk1 || min_dist_k > nk2 )
  {
    is_here = 0;
    *si = -1000;
    *sk = -1000;
    *sx_inc = 0.0;
    *sz_inc = 0.0;
    return is_here;
  }

  // in this thread
  is_here = 1;

  float points_x[4];
  float points_z[4];
  float points_i[4];
  float points_k[4];

  for (int kk=0; kk<2; kk++)
  {
      for (int ii=0; ii<2; ii++)
      {
        int cur_i = min_dist_i-1+ii;
        int cur_k = min_dist_k-1+kk;

        for (int n3=0; n3<2; n3++) {
            for (int n1=0; n1<2; n1++) {
              int iptr_cube = n1 + n3 * 2;
              int iptr = (cur_i+n1) + (cur_k+n3) * siz_iz;
              points_x[iptr_cube] = x2d[iptr];
              points_z[iptr_cube] = z2d[iptr];
              points_i[iptr_cube] = cur_i+n1;
              points_k[iptr_cube] = cur_k+n3;
            }
        }

        if (fdlib_math_isPoint2InQuad(sx,sz,points_x,points_z) == 1)
        {
          float si_curv, sk_curv;

          gd_curv_coord2index_sample(sx,sz,
              points_x,points_z,
              points_i,points_k,
              100,100,
              &si_curv, &sk_curv);

          // convert to return values
          *si = min_dist_i;
          *sk = min_dist_k;
          *sx_inc = si_curv - min_dist_i;
          *sz_inc = sk_curv - min_dist_k;

          return is_here;
        }
      }
  }

  // if not in any cube due to bug, set default value
  //    if everything is right, should be return 10 line before
  *si = -1000;
  *sk = -1000;
  *sx_inc = 0.0;
  *sz_inc = 0.0;

  return is_here;
}

/*
 * convert depth to axis
 */
int
gd_curv_depth_to_axis(gd_t *gd,
                      float sx,
                      float *sz)
{
  int ierr = 0;
  int nk2 = gd->nk2;
  int ni1 = gd->ni1;
  int ni2 = gd->ni2;
  size_t siz_iz = gd->siz_iz;
  size_t iptr;
  size_t iptr_pt;
  float points_x[2];
  float points_z[2];
  float xmin, xmax;

  // not here if outside coord range
  if ( sx < gd->xmin || sx > gd->xmax )
  {
    return ierr;
  }


  for(int i=ni1; i<ni2; i++)
  {
    iptr = i + nk2 * siz_iz;
    xmin = gd->x2d[iptr];
    xmax = gd->x2d[iptr+1];
    if(sx>=xmin  && sx <xmax )
    {
      for (int n1=0; n1<2; n1++) 
      {
        iptr_pt = (i+n1) + nk2 * siz_iz;
        points_x[n1] = gd->x2d[iptr_pt];
        points_z[n1] = gd->z2d[iptr_pt];
      }
    }
  }

  float ztopo = linear_interp_1d(sx,points_x,points_z);
  *sz = ztopo - (*sz);

  return ierr;
}

float
linear_interp_1d(float ix, float *x, float *z)
{
  float iz;
  iz = (ix-x[1])/(x[0]-x[1])*z[0] + (ix-x[0])/(x[1]-x[0]) * z[1];

  return iz;
}


/* 
 * find curv index using sampling
 */

int
gd_curv_coord2index_sample(float sx, float sz, 
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    int    nx_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sk_curv)
{
  float Lx[2], Lz[2];

  // init closest point
  float min_dist = sqrtf(  (sx - points_x[0]) * (sx - points_x[0])
      + (sz - points_z[0]) * (sz - points_z[0]) );
  int min_dist_i = 0 ;
  int min_dist_k = 0 ;

  // linear interp for all sample
  for (int n3=0; n3<nz_sample+1; n3++)
  {
    Lz[1] = (float)(n3) / (float)(nz_sample);
    Lz[0] = 1.0 - Lz[1];
      for (int n1=0; n1<nx_sample+1; n1++)
      {
        Lx[1] = (float)(n1) / (float)(nx_sample);
        Lx[0] = 1.0 - Lx[1];

        // interp
        float x_pt=0;
        float z_pt=0;
        for (int kk=0; kk<2; kk++) {
            for (int ii=0; ii<2; ii++)
            {
              int iptr_cube = ii + kk * 2;
              x_pt += Lx[ii]*Lz[kk] * points_x[iptr_cube];
              z_pt += Lx[ii]*Lz[kk] * points_z[iptr_cube];
            }
        }

        // find min dist
        float DistInt = sqrtf(  (sx - x_pt) * (sx - x_pt)
            + (sz - z_pt) * (sz - z_pt) );

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = n1;
          min_dist_k = n3;
        }
      } // n1
  } // n3

  *si_curv = points_i[0] + (float)min_dist_i / (float)nx_sample;
  *sk_curv = points_k[0] + (float)min_dist_k / (float)nz_sample;

  return 0;
}

/* 
 * interp curv coord using inverse distance interp
 */

int
gd_curv_coord2index_rdinterp(float sx, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    float *si_curv, // interped curv coord
    float *sk_curv)
{
  float weight[num_points];
  float total_weight = 0.0 ;

  // cal weight
  int at_point_indx = -1;
  for (int i=0; i<num_points; i++)
  {
    float dist = sqrtf ((sx - points_x[i]) * (sx - points_x[i])
        + (sz - points_z[i]) * (sz - points_z[i])
        );
    if (dist < 1e-9) {
      at_point_indx = i;
    } else {
      weight[i]   = 1.0 / dist;
      total_weight += weight[i];
    }
  }
  // if at a point
  if (at_point_indx > 0) {
    total_weight = 1.0;
    // other weight 0
    for (int i=0; i<num_points; i++) {
      weight[i] = 0.0;
    }
    // point weight 1
    weight[at_point_indx] = 1.0;
  }

  // interp

  *si_curv = 0.0;
  *sk_curv = 0.0;

  for (int i=0; i<num_points; i++)
  {
    weight[i] *= 1.0 / total_weight ;

    (*si_curv) += weight[i] * points_i[i];
    (*sk_curv) += weight[i] * points_k[i];  

    fprintf(stdout,"---- i=%d,weight=%f,points_i=%f,points_k=%f\n",
        i,weight[i],points_i[i],points_k[i]);
  }

  return 0;
}

float
gd_coord_get_x(gd_t *gd, int i, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->x1d[i];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + k * gd->siz_iz;
    var = gd->x2d[iptr];
  }

  return var;
}

float
gd_coord_get_z(gd_t *gd, int i, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->z1d[k];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + k * gd->siz_iz;
    var = gd->z2d[iptr];
  }

  return var;
}

//
// set grid size
//

int
gd_indx_set(gd_t *const gd,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_z,
            const int fdx_nghosts,
            const int fdz_nghosts,
            const int verbose)
{
  int ierr = 0;

  int ni = number_of_total_grid_points_x;
  int nk = number_of_total_grid_points_z;
  
  // add ghost points
  int nx = ni + 2 * fdx_nghosts;
  int nz = nk + 2 * fdz_nghosts;

  gd->ni = ni;
  gd->nk = nk;

  gd->nx = nx;
  gd->nz = nz;

  gd->ni1 = fdx_nghosts;
  gd->ni2 = gd->ni1 + ni - 1;

  gd->nk1 = fdz_nghosts;
  gd->nk2 = gd->nk1 + nk - 1;

  gd->gni1 = 0;
  gd->gni2 = gd->gni1 + ni - 1;

  gd->gnk1 = 0;
  gd->gnk2 = gd->gnk1 + nk - 1;

  gd->siz_iz  = gd->nx;
  gd->siz_icmp = gd->nx * gd->nz;

  // set npoint_ghosts according to fdz_nghosts
  gd->npoint_ghosts = fdz_nghosts;

  gd->fdx_nghosts = fdx_nghosts;
  gd->fdz_nghosts = fdz_nghosts;

  gd->index_name = fdlib_mem_malloc_2l_char(
                        CONST_NDIM, CONST_MAX_STRLEN, "gd name");

  // grid coord name
  sprintf(gd->index_name[0],"%s","i");
  sprintf(gd->index_name[1],"%s","k");

  return ierr;
}

/*
 * give a local index ref, check if in this thread
 */

int
gd_lindx_is_inner(int i, int k, gd_t *gd)
{
  int is_in = 0;

  if (   i >= gd->ni1 && i <= gd->ni2
      && k >= gd->nk1 && k <= gd->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_pindx_is_inner(int i_phy, int k_phy, gd_t *gd)
{
  int is_in = 0;

  int i = i_phy + gd->fdx_nghosts;
  int k = k_phy + gd->fdz_nghosts;

  if (   i >= gd->ni1 && i <= gd->ni2
      && k >= gd->nk1 && k <= gd->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_pindx_is_inner_i(int i_phy, gd_t *gd)
{
  int is_in = 0;

  int i = i_phy + gd->fdx_nghosts;

  if (   i >= gd->ni1 && i <= gd->ni2)
  {
    is_in = 1;
  }

  return is_in;
}  

int
gd_pindx_is_inner_k(int k_phy, gd_t *gd)
{
  int is_in = 0;

  int k = k_phy + gd->fdz_nghosts;

  if ( k >= gd->nk1 && k <= gd->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

/*
 * print for QC
 */

int
gd_print(gd_t *gd)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " nx    = %-10d\n", gd->nx);
  fprintf(stdout, " nz    = %-10d\n", gd->nz);
  fprintf(stdout, " ni    = %-10d\n", gd->ni);
  fprintf(stdout, " nk    = %-10d\n", gd->nk);

  fprintf(stdout, " ni1   = %-10d\n", gd->ni1);
  fprintf(stdout, " ni2   = %-10d\n", gd->ni2);
  fprintf(stdout, " nk1   = %-10d\n", gd->nk1);
  fprintf(stdout, " nk2   = %-10d\n", gd->nk2);

  return 0;
}



