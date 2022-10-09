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
   * 0-2: x3d, y3d, z3d
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
    
  // extend to ghosts. may replaced by mpi exchange
  // x1, mirror
  for (size_t k = 0; k < nz; k++){
      for (size_t i = 0; i < ni1; i++)
      {
        size_t iptr = i + k * siz_iz;
        jac2d[iptr] = jac2d[iptr + (ni1-i)*2 -1 ];
         xi_x[iptr] =  xi_x[iptr + (ni1-i)*2 -1 ];
         xi_z[iptr] =  xi_z[iptr + (ni1-i)*2 -1 ];
         zt_x[iptr] =  zt_x[iptr + (ni1-i)*2 -1 ];
         zt_z[iptr] =  zt_z[iptr + (ni1-i)*2 -1 ];
      }
  }
  // x2, mirror
  for (size_t k = 0; k < nz; k++){
      for (size_t i = ni2+1; i < nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        jac2d[iptr] = jac2d[iptr - (i-ni2)*2 +1 ];
         xi_x[iptr] =  xi_x[iptr - (i-ni2)*2 +1 ];
         xi_z[iptr] =  xi_z[iptr - (i-ni2)*2 +1 ];
         zt_x[iptr] =  zt_x[iptr - (i-ni2)*2 +1 ];
         zt_z[iptr] =  zt_z[iptr - (i-ni2)*2 +1 ];
      }
  }
  // z1, mirror
  for (size_t k = 0; k < nk1; k++) {
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + k * siz_iz;
        jac2d[iptr] = jac2d[iptr + ((nk1-k)*2 -1) * siz_iz ];
         xi_x[iptr] =  xi_x[iptr + ((nk1-k)*2 -1) * siz_iz ];
         xi_z[iptr] =  xi_z[iptr + ((nk1-k)*2 -1) * siz_iz ];
         zt_x[iptr] =  zt_x[iptr + ((nk1-k)*2 -1) * siz_iz ];
         zt_z[iptr] =  zt_z[iptr + ((nk1-k)*2 -1) * siz_iz ];
      }
  }
  // z2, mirror
  for (size_t k = nk2+1; k < nz; k++) {
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + k * siz_iz;
        jac2d[iptr] = jac2d[iptr - ((k-nk2)*2 -1) * siz_iz ];
         xi_x[iptr] =  xi_x[iptr - ((k-nk2)*2 -1) * siz_iz ];
         xi_z[iptr] =  xi_z[iptr - ((k-nk2)*2 -1) * siz_iz ];
         zt_x[iptr] =  zt_x[iptr - ((k-nk2)*2 -1) * siz_iz ];
         zt_z[iptr] =  zt_z[iptr - ((k-nk2)*2 -1) * siz_iz ];
      }
  }

  return;
}

/*
 * generate cartesian grid for curv struct
 */
void
gd_curv_gen_cart(gd_t *gdcurv,
                 float dx, float x0_glob,
                 float dz, float z0_glob)
{
  float *x3d = gdcurv->x2d;
  float *z3d = gdcurv->z2d;

  float x0 = x0_glob + (0 - gdcurv->fdx_nghosts) * dx;
  float z0 = z0_glob + (0 - gdcurv->fdz_nghosts) * dz;

  size_t iptr = 0;
  for (size_t k=0; k<gdcurv->nz; k++)
  {
      for (size_t i=0; i<gdcurv->nx; i++)
      {
        x3d[iptr] = x0 + i * dx;
        z3d[iptr] = z0 + k * dz;

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
   * 0-2: x3d, y3d, z3d
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
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
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

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v3d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

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

  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"open coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // read vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++)
  {
    float *ptr = gdcurv->v3d + gdcurv->cmp_pos[ivar];

    ierr = nc_inq_varid(ncid, gdcurv->cmp_name[ivar], &varid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }

    ierr = nc_get_var(ncid, varid, ptr);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}

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
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, 1, dimid+1, &varid[0]);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, 1, dimid+0, &varid[1]);

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int l_count[] = { ni, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  ierr = nc_put_var_float(ncid, varid[0], gdcart->x1d);
  ierr = nc_put_var_float(ncid, varid[1], gdcart->z1d);
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

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
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, g3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
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

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = metric->v3d + g3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
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
  if (ierr != NC_NOERR){
    fprintf(stderr,"open coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // read vars
  for (int ivar=0; ivar<metric->ncmp; ivar++)
  {
    float *ptr = metric->v3d + metric->cmp_pos[ivar];

    ierr = nc_inq_varid(ncid, metric->cmp_name[ivar], &varid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }

    ierr = nc_get_var(ncid, varid, ptr);
    if (ierr != NC_NOERR){
      fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}

/*
 * generate grid from a input grid layer file
 * code by SunWenliang  2021/06/14
 */
int
gd_curv_gen_layer(char *in_grid_layer_file,
                      int *grid_layer_resample_factor,
                      int *grid_layer_start,
                      int n_total_grid_x,
                      int n_total_grid_z,
                      float * x2d,
                      float * z2d,
                      int nx, int ni, int gni1, int fdx_nghosts, 
                      int nz, int nk, int gnk1, int fdz_nghosts)
{
  int ierr = 0;

  int nLayers;
  int n_Interfaces;
  int nx_layers;
  int nx_first;
  int nz_first;
  int nx_interp;
  int nz_interp;
  size_t iptr1 ;
  size_t iptr2 ;

//
// read data from layer file
//

  FILE * fp;
  fp = fopen( in_grid_layer_file, "r");
  if (!fp){
    fprintf(stderr,"Failed to open input file of interface!!\n");
    exit(-1);
  }

  // read number of interface
  fscanf(fp,"%d", &n_Interfaces);

  nLayers = n_Interfaces - 1;
  int NCellPerlay_IN[nLayers];
  int NCellPerlay[nLayers];
  int VmapSpacingIsequal[nLayers];

  // read cells and is_equal of each layer
  for ( int i=0; i<nLayers; i++)
  {
    fscanf(fp,"%d",&NCellPerlay_IN[i]);
  }
  for ( int i=0; i<nLayers; i++)
  {
    fscanf(fp,"%d",&VmapSpacingIsequal[i]);
  }

  // read number of horizontal sampling
  fscanf(fp,"%d",&nx_layers);

  size_t siz_icmp_layerIn = nx_layers * (nLayers+1) ;  
  float * layer2d_In        = NULL;
  layer2d_In = (float *)fdlib_mem_malloc_1d(siz_icmp_layerIn * 2 * sizeof(float), 
                                                              "gd_curv_gen_layer");
  for ( int i=0; i<siz_icmp_layerIn; i++)
  {
    // read x,z of each sample
    fscanf(fp,"%f",&layer2d_In[i                     ]);
    fscanf(fp,"%f",&layer2d_In[i + siz_icmp_layerIn ]);
  }
  fclose( fp );

//
// resample of input interface grid nodes
//

  int x_interp_factor = abs(grid_layer_resample_factor[0]);
  int z_interp_factor = abs(grid_layer_resample_factor[1]);

  // effective layer nx after downsampling
  if ( grid_layer_resample_factor[0] < 0 ) {
    nx_interp = floor((nx_layers              + 1) / x_interp_factor);
    nx_first  = floor((grid_layer_start[0] + 1) / x_interp_factor);
  }
  // effective layer nx after upsampling
  else
  { 
    nx_interp = (nx_layers               -1) * x_interp_factor + 1;
    nx_first  = (grid_layer_start[0]-1) * x_interp_factor + 1;
  }

  nz_interp = 1;
  if ( grid_layer_resample_factor[1] < 0 ) 
  {
    for ( int kk=0; kk<nLayers; kk++)
    {
      NCellPerlay[kk] = NCellPerlay_IN[kk] / z_interp_factor ;
      nz_interp = nz_interp + NCellPerlay[kk];
    }
  }
  else
  { 
    for ( int kk=0; kk<nLayers; kk++)
    {
      NCellPerlay[kk] = NCellPerlay_IN[kk] * z_interp_factor ;
      nz_interp = nz_interp + NCellPerlay[kk];
    }
  }

  nz_first = nz_interp - grid_layer_start[1] - n_total_grid_z;
  int x_gd_first = nx_first + gni1 - fdx_nghosts; //Starting x index at input interface
  int z_gd_first = nz_first + gnk1 - fdz_nghosts;

  if (gni1== 0)
  {
    fprintf(stdout, "-------------------------------------------------------\n");
    fprintf(stdout, "--> Input layer information.\n");
    fprintf(stdout, "-------------------------------------------------------\n");
    fprintf(stdout, "n_Interfaces = %d\n", n_Interfaces);
    fprintf(stdout, "nx_layers = %d\n", nx_layers);

    fprintf(stdout, "     resample_factor = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", grid_layer_resample_factor[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "         NCellPerlay = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", NCellPerlay_IN[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "NCellPerlay_resample = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", NCellPerlay[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "  VmapSpacingIsequal = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", VmapSpacingIsequal[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "n_input_resample = [ %d %d ]\n", nx_interp, nz_interp);
    fprintf(stdout, "          n_c3d  = [ %d %d ]\n", n_total_grid_x + fdx_nghosts * 2, 
                                                      n_total_grid_z + fdz_nghosts * 2);
    fprintf(stdout, "horizontal_index = [ %d ] --> [ %d ]\n", x_gd_first,  
                                               x_gd_first + n_total_grid_x + fdx_nghosts * 2 );
    fprintf(stdout, "  vertical_index = %d --> %d \n", z_gd_first, z_gd_first + n_total_grid_z+ 
                                                                           fdz_nghosts);
    fprintf(stdout, "-------------------------------------------------------\n");
    fflush(stdout);
  }

  //chech interface parameters
  if (nx_first - fdx_nghosts < 0)
  {
    fprintf(stdout, "Input Parameter Error: horizontal_start_index.\n");
    fflush(stdout);
    fprintf(stdout, "horizontal_start_index is related to fdx_nghosts\n");
    fflush(stdout);
    exit(-1);
  }
  else if (nx_first + n_total_grid_x + fdx_nghosts > nx_interp)
  {
    fprintf(stdout, "Input Parameter Error:\n"); 
    fflush(stdout);
    fprintf(stdout, "please check the number_of_total_grid_points， refine_factor and in_grid_layer_file file !!\n"); 
    fflush(stdout);
    exit(-1);
  }
  else if ( z_gd_first < 0)
  {
    fprintf(stdout, "Input Parameter Error:\n"); 
    fflush(stdout);
    fprintf(stdout, "please check the vertical_FreeSurf_AfterRefine_start_index\n number_of_total_grid_point_z\n refine_factor \nin_grid_layer_file file\n"); 
    fflush(stdout);
    exit(-1);
  }

  // resample based on interp_factor on input grid layer
  size_t siz_icmp_layer_interp = nx_interp * (nLayers + 1);
  float *layer2d_interp = NULL;
  layer2d_interp = (float *)fdlib_mem_malloc_1d(siz_icmp_layer_interp * 2 * sizeof(float),"gd_curv_gen_layer"); 
  
  // 
  float   *xi_lenX = NULL;
  float   *xn_lenX = NULL;
  float *yi_x_lenX = NULL;
  float *yn_x_lenX = NULL;
  float *yi_z_lenX = NULL;
  float *yn_z_lenX = NULL;

  xi_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  xn_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");
  yi_x_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  yn_x_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");
  yi_z_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  yn_z_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");

  // Downsampling in X 
  if (grid_layer_resample_factor[0] < 2) 
  {
    for (int kk = 0; kk < nLayers+1; kk++) 
    {
      for (int ii = 0; ii < nx_interp; ii++)
      {
        iptr1 = M_gd_INDEX(ii, kk, nx_interp );
        iptr2 = M_gd_INDEX(ii * x_interp_factor,  kk, nx_layers);

        layer2d_interp[ iptr1                        ] = layer2d_In[ iptr2 ];
        layer2d_interp[ iptr1 +siz_icmp_layer_interp] = layer2d_In[ iptr2 +siz_icmp_layerIn]; 
      }
    }
  }
  //Interpolating in in X
  else if( grid_layer_resample_factor[0] >=2 )  
  {
    for ( int ii1=0; ii1<nx_interp; ii1++) xi_lenX[ii1] = ii1;
    for ( int ii1=0; ii1<nx_layers; ii1++) xn_lenX[ii1] = ii1*x_interp_factor;

    for (int kk = 0; kk < nLayers+1; kk++) {
      for ( int ii1=0; ii1<nx_layers; ii1++)
      {
        yn_x_lenX[ii1] = layer2d_In[M_gd_INDEX(ii1, kk, nx_layers)                     ];
        yn_z_lenX[ii1] = layer2d_In[M_gd_INDEX(ii1, kk, nx_layers) + siz_icmp_layerIn ];
      }

      gd_SPL(nx_layers, xn_lenX, yn_x_lenX, nx_interp, xi_lenX, yi_x_lenX);
      gd_SPL(nx_layers, xn_lenX, yn_z_lenX, nx_interp, xi_lenX, yi_z_lenX);

      for (int ii = 0; ii < nx_interp; ii++)
      {
          iptr1 = M_gd_INDEX(ii, kk, nx_interp);
          layer2d_interp[iptr1] = yi_x_lenX[ii];
          layer2d_interp[iptr1 + siz_icmp_layer_interp] = yi_z_lenX[ii];
      }
    }
  }

  free(xi_lenX);
  free(xn_lenX);
  free(yi_x_lenX);
  free(yn_x_lenX);
  free(yi_z_lenX);
  free(yn_z_lenX);

//
// generate FD grid
//

  size_t siz_slcie = nx * nz;
  float *layer2d = NULL;
  layer2d = ( float * ) malloc( sizeof( float ) * nx * (nLayers+1) * 2 );

  // suffix ch means:
  float zdiff_two_interface; // thickness between interfaces
  float zdiff_two_interface_ch;

  float *xlayer2d = layer2d + 0 * nx*(nLayers+1);
  float *zlayer2d = layer2d + 1 * nx*(nLayers+1);
  
  float * zlayerpart  = NULL;
  float * xlayerpart  = NULL;
  zlayerpart = (float *)fdlib_mem_malloc_1d((nLayers+1) * sizeof(float), "gd_curv_gen_layer");
  xlayerpart = (float *)fdlib_mem_malloc_1d((nLayers+1) * sizeof(float), "gd_curv_gen_layer");

  float * z2dpart  = NULL;
  float * x2dpart  = NULL;
  z2dpart = (float *)fdlib_mem_malloc_1d((nz_interp+fdz_nghosts) * sizeof(float), "gd_curv_gen_layer");
  x2dpart = (float *)fdlib_mem_malloc_1d((nz_interp+fdz_nghosts) * sizeof(float), "gd_curv_gen_layer");

  // layer regions of input model covered by this thread
  for (int k = 0; k < (nLayers + 1); k++) {
    for (int i = 0; i < nx; i++)
    {
      iptr1 = M_gd_INDEX(i, k, nx);
      iptr2 = M_gd_INDEX(i + x_gd_first, k, nx_interp);
      xlayer2d[iptr1] = layer2d_interp[iptr2];
      zlayer2d[iptr1] = layer2d_interp[iptr2 + siz_icmp_layer_interp];
    }
  }

  // interp both z and x,y
  for (int i = 0; i < nx; i++)
  {
    for (int k = 0; k < nLayers + 1; k++)
    {
      xlayerpart[k] = xlayer2d[M_gd_INDEX(i, k, nx)];
      zlayerpart[k] = zlayer2d[M_gd_INDEX(i, k, nx)];
    }
    // Interpolating the Z coordinate of the grid 
    gd_grid_z_interp(z2dpart, zlayerpart, NCellPerlay, VmapSpacingIsequal,
                     nLayers);

    // Interpolating the X and Y coordinate of the grid ...
    // by cubic spline interpolation method.
    gd_SPL(nLayers + 1, zlayerpart, xlayerpart, nz_interp, z2dpart, x2dpart);

    //Grids outside the Free Surface.
    for (int k = 1; k < fdz_nghosts + 1; k++)
    {
      iptr1 = nz_interp - 1;
      x2dpart[iptr1 + k] = x2dpart[iptr1] * 2 - x2dpart[iptr1 - k];
      z2dpart[iptr1 + k] = z2dpart[iptr1] * 2 - z2dpart[iptr1 - k];
    }

    for ( int k = 0; k < nk + fdz_nghosts*2; k ++)
    {
      iptr1 = z_gd_first + k;
      x2d[M_gd_INDEX(i, k, nx)] = x2dpart[iptr1];
      z2d[M_gd_INDEX(i, k, nx)] = z2dpart[iptr1];
    }
  }

  free(layer2d);
  free(layer2d_In);
  free(layer2d_interp);
  free(zlayerpart);
  free(xlayerpart);
  free(z2dpart);
  free(x2dpart);

  return ierr;
}

//  Interpolating the Z coordinate
int gd_grid_z_interp(float *z3dpart, float *zlayerpart, int *NCellPerlay,
                     int *VmapSpacingIsequal, int nLayers)
{
  int ierr = 0;

  // If Grid step size of Z coordinate is not equal, divide the layer into two parts.
  // N1 is first part number of this layer and N2 is second.
  int N1;
  int N2;
  float distance;  //Interface spacing
  float dmidpoint;  //The length of the N1th grid
  float N1_Isochromatic; 
  float N2_Isochromatic;
  int NCellPerlay_max = 0;
  int sumNCellPerlay = 0;
  float LayerDz[nLayers];

  for (int i = 0; i < nLayers; i++)
  {
    if (NCellPerlay_max < NCellPerlay[i])
    {
      NCellPerlay_max = NCellPerlay[i];
    }
  }

  float range1[NCellPerlay_max];
  float range2[NCellPerlay_max];

  for (int i = 0; i < nLayers; i++)
  {
    LayerDz[i] = (zlayerpart[i + 1] - zlayerpart[i]) / NCellPerlay[i];
  }

  for (int i = 0; i < nLayers; i++)
  {
    // The grid spacing is equal
    z3dpart[sumNCellPerlay] = zlayerpart[i];
    for (int ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i] - 1; ii++)
    {
      z3dpart[ii+1] = zlayerpart[i] + (ii - sumNCellPerlay + 1) * LayerDz[i];
    }
    // The grid spacing is not equal
    if (VmapSpacingIsequal[i] < 1 && i > 0 && i < nLayers - 1)
    {
      range1[0] = zlayerpart[i];
      if ((LayerDz[i] - LayerDz[i - 1]) * (LayerDz[i + 1] - LayerDz[i]) > 0)
      {
        N1 = ((LayerDz[i] - LayerDz[i - 1]) / (LayerDz[i + 1] - LayerDz[i - 1])) * (NCellPerlay[i] + 2) + 1;
        if (N1 > NCellPerlay[i] - 1)
          N1 = NCellPerlay[i] - 1;
        N2 = NCellPerlay[i] + 2 - N1;

        distance = zlayerpart[i + 1] - zlayerpart[i] + LayerDz[i - 1] + LayerDz[i + 1];
        dmidpoint = (2 * distance - LayerDz[i - 1] * N1 - LayerDz[i + 1] * (N2 + 1)) / (N1 + N2 - 1);

        N1_Isochromatic = (dmidpoint - LayerDz[i - 1]) / (N1 - 1);
        N2_Isochromatic = (LayerDz[i + 1] - dmidpoint) / (N2);

        for (int j = 1; j < N1; j++)
        {
          range1[j] = range1[j - 1] + LayerDz[i - 1] + j * N1_Isochromatic;
        }
        for (int j = N1; j < N1 + N2; j++)
        {
          range1[j] = range1[j - 1] + LayerDz[i + 1] - (N2 + N1 - j - 1) * N2_Isochromatic;
        }
        for (int ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii++)
        {
          z3dpart[ii] = range1[ii - sumNCellPerlay];
        }
      }
    }
    sumNCellPerlay += NCellPerlay[i];
  }
  z3dpart[sumNCellPerlay] = zlayerpart[nLayers];

  return ierr;
}

float gd_seval(int ni, float u,
            int n, float x[], float y[],
            float b[], float c[], float d[],
            int *last)
{
  int i, j, k;
  float w;
  i = *last;
  if (i >= n - 1) i = 0;
  if (i < 0) i = 0;
  if ((x[i] > u) || (x[i + 1] < u))//??
  {
    i = 0;
    j = n;
    do
    {
      k = (i + j) / 2;
      if (u < x[k]) j = k;
      if (u >= x[k]) i = k;
    } while (j > i + 1);
  }
  *last = i;
  w = u - x[i];
  w = y[i] + w * (b[i] + w * (c[i] + w * d[i]));
  return (w);
}

//// Cubic spline difference function
int gd_SPLine( int n, int end1, int end2,
           float slope1, float slope2,
           float x[], float y[],
           float b[], float c[], float d[],
           int *iflag)
{
  int nm1, ib, i, ascend;
  float t;
  nm1 = n - 1;
  *iflag = 0;
  if (n < 2)
  {
    *iflag = 1;
    goto Leavegd_SPLine;
  }
  ascend = 1;
  for (i = 1; i < n; ++i) if (x[i] <= x[i - 1]) ascend = 0;
  if (!ascend)
  {
    *iflag = 2;
    goto Leavegd_SPLine;
  }
  if (n >= 3)
  {
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (i = 1; i < nm1; ++i)
    {
      d[i] = x[i + 1] - x[i];
      b[i] = 2.0 * (d[i - 1] + d[i]);
      c[i + 1] = (y[i + 1] - y[i]) / d[i];
      c[i] = c[i + 1] - c[i];
    }
    b[0] = -d[0];
    b[nm1] = -d[n - 2];
    c[0] = 0.0;
    c[nm1] = 0.0;
    if (n != 3)
    {
      c[0] = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
      c[nm1] = c[n - 2] / (x[nm1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4]);
      c[0] = c[0] * d[0] * d[0] / (x[3] - x[0]);
      c[nm1] = -c[nm1] * d[n - 2] * d[n - 2] / (x[nm1] - x[n - 4]);
    }
    if (end1 == 1)
    {
      b[0] = 2.0 * (x[1] - x[0]);
      c[0] = (y[1] - y[0]) / (x[1] - x[0]) - slope1;
    }
    if (end2 == 1)
    {
      b[nm1] = 2.0 * (x[nm1] - x[n - 2]);
      c[nm1] = slope2 - (y[nm1] - y[n - 2]) / (x[nm1] - x[n - 2]);
    }
    // Forward elimination 
    for (i = 1; i < n; ++i)
    {
      t = d[i - 1] / b[i - 1];
      b[i] = b[i] - t * d[i - 1];
      c[i] = c[i] - t * c[i - 1];
    }
    // Back substitution 
    c[nm1] = c[nm1] / b[nm1];
    for (ib = 0; ib < nm1; ++ib)
    {
      i = n - ib - 2;
      c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
    }
    b[nm1] = (y[nm1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[nm1]);
    for (i = 0; i < nm1; ++i)
    {
      b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
      d[i] = (c[i + 1] - c[i]) / d[i];
      c[i] = 3.0 * c[i];
    }
    c[nm1] = 3.0 * c[nm1];
    d[nm1] = d[n - 2];
  }
  else
  {
    b[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c[0] = 0.0;
    d[0] = 0.0;
    b[1] = b[0];
    c[1] = 0.0;
    d[1] = 0.0;
  }
Leavegd_SPLine:
  return 0;
}

//// Cubic spline difference function in grid interpolation
void gd_SPL(int n, float *x, float *y, int ni, float *xi, float *yi)
{
  float *b, *c, *d;
  int iflag=0, last=0, i=0;
  b = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");
  c = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");
  d = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");

  if (!d) { printf("no enough memory for b,c,d\n"); }
  else {
    gd_SPLine(n, 0, 0, 0.0, 0.0, x, y, b, c, d, &iflag);
    for (i = 0; i<ni; i++)
      yi[i] = gd_seval(ni, xi[i], n, x, y, b, c, d, &last);
    free(b);
    free(c);
    free(d);
  };
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

  float * x3d = gd->x2d;
  float * z3d = gd->z2d;

  // outside coord range
  if ( sx < gd->xmin || sx > gd->xmax ||
       sz < gd->zmin || sz > gd->zmax)
  {
    is_here = 0;
    return is_here;
  }

  // init closest point
  float min_dist = sqrtf(  (sx - x3d[0]) * (sx - x3d[0])
      + (sz - z3d[0]) * (sz - z3d[0]) );
  int min_dist_i = 0 ;
  int min_dist_k = 0 ;

  // compute distance to each point
  for (int k=0; k<nz; k++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float x = x3d[iptr];
        float z = z3d[iptr];

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
              int iptr = (cur_i+n1)  + (cur_k+n3) * siz_iz;
              points_x[iptr_cube] = x3d[iptr];
              points_z[iptr_cube] = z3d[iptr];
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
  *si = min_dist_i;
  *sk = min_dist_k;
  *sx_inc = 0.0;
  *sz_inc = 0.0;

  return is_here;
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



