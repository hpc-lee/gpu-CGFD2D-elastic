/******************************************************************************
 *
 * This function is used to discretize a given grid model (binary version) 
 *  to a calculation grid.
 *
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    03/2022: Created by Luqian Jiang 
 *
 *******************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <math.h>
#include "media_geometry2d.hpp"
#include "media_bin2model.hpp"
#include "media_utility.hpp"
#include "media_read_file.hpp"


/* =================== for C call======================== */
int media_bin2model_el_iso(
    float *rho2d,
    float *lam2d,
    float *mu2d, 
    const float *x2d,
    const float *z2d,
    size_t nx, size_t nz,
    float xmin, float xmax,
    int grid_type,
    int  *bin_order,    // eg, [1, 0]=[z, x] 0:x, 1:z
    int  *bin_size,     // [ndim1, ndim2],
    float  *bin_spacing,  // [dh1, dh2],
    float  *bin_origin,   // [h0_1, h0_2],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs  )
{
  
    // get the dim order and set error reporting
    int dimx = -1, dimz = -1;

    for (int i = 0; i < 2; i++) {
      if (bin_order[i] == 0) dimx = i;
      else if (bin_order[i] == 1) dimz = i;
      else {
        fprintf(stderr, "Error: wrong number in bin_order: bin_order[%d] = %d! \n"\
                        "       %d is a wrong order, it must be 0, 1, or 2. \n",
                i, bin_order[i], bin_order[i] );
        fflush(stderr);
        exit(1);
      }
    }

    if (dimx == -1) {
      fprintf(stderr, "Error: missing the order of x in bin_order.\n");
      fprintf(stderr, "       bin_order=%d,%d,%d\n", bin_order[0],bin_order[1]);
      fflush(stderr);
      exit(1);
    }
    if (dimz == -1) {
      fprintf(stderr, "Error: missing the order of z in bin_order.\n");
      fprintf(stderr, "       bin_order=%d,%d,%d\n", bin_order[0],bin_order[1]);
      fflush(stderr);
      exit(1);
    }
    
    // get the read range index
    // x should be increasing
    float  bin_x0 = bin_origin[dimx]; 
    float  bin_z0 = bin_origin[dimz]; 
    float  bin_dx = bin_spacing[dimx]; 
    float  bin_dz = bin_spacing[dimz];

    int ix_start = (xmin < bin_x0 ? 0:( xmin - bin_x0 ) / bin_dx); 
    if (ix_start >= bin_size[dimx] || xmax < bin_x0) {
      fprintf(stderr, "Error: The given model does not coincide with the calculation area\n"\
                      "       in x-direction. Please check the model settings!");
      fflush(stderr);
      exit(1);
    }
    int ix_end = ceil( ( xmax - bin_x0 ) / bin_dx );
    if (ix_end >= bin_size[dimx])
      ix_end = bin_size[dimx]-1;

    // judge whether the calculation area does coincide with the model in z direction
    float zmin = FLT_MAX, zmax = -FLT_MAX;
    size_t siz_slice =  nx * nz;
    size_t siz_z2d = (grid_type == GRID_CART ? nz:siz_slice);

    for (size_t i = 0; i < siz_z2d; i++) {
      zmin = std::min(zmin, z2d[i]);
      zmax = std::max(zmax, z2d[i]);
    }
    
    float bin_zn = bin_z0 + bin_dz * (bin_size[dimz]-1);
    if (zmin >= std::max(bin_zn, bin_z0) || zmax <= std::min(bin_z0, bin_zn) ) {
      fprintf(stderr, "Error: The given model does not coincide with the calculation area\n"\
                      "       in z-direction. Please check the model settings!");
      fflush(stderr);
      exit(1);
    }

    // the coordinate vectors of the given model (after cutting by [min max] )
    size_t bin_nx = ix_end-ix_start+1;
    size_t bin_nz = bin_size[dimz];
    std::vector<float> xvec(bin_nx);
    std::vector<float> zvec(bin_nz);

    for (size_t ix = 0; ix < bin_nx; ix++) {
      xvec[ix] = bin_x0 + (ix_start+ix)*bin_dx;
    }

    for (size_t iz = 0; iz < bin_nz; iz++) {
      zvec[iz] = bin_z0 + iz*bin_dz;
    }

    int bin_start[2];
    int bin_end[2];
    
    bin_start[dimx] = ix_start; bin_end[dimx] = ix_end;
    bin_start[dimz] = 0; bin_end[dimz] = bin_nz-1;

    //- Read bin file 
    size_t bin_line = bin_nx;
    size_t bin_slice = bin_nx * bin_nz;
    float *bin_rho = new float[bin_slice];
    float *bin_vp  = new float[bin_slice];
    float *bin_vs  = new float[bin_slice];

    fprintf(stdout, "- reading model file: %s, \n", bin_file_rho);
    read_bin_file(bin_file_rho, bin_rho, dimx, dimz, bin_start, bin_end, bin_size, bin_line);

    fprintf(stdout, "                      %s, \n", bin_file_vp);
    read_bin_file(bin_file_vp,  bin_vp,  dimx, dimz, bin_start, bin_end, bin_size, bin_line);

    fprintf(stdout, "                      %s, \n", bin_file_vs);
    read_bin_file(bin_file_vs,  bin_vs,  dimx, dimz, bin_start, bin_end, bin_size, bin_line);

    // media parameterization
    parameterization_bin_el_iso_loc(rho2d, lam2d, mu2d, x2d, z2d, nx, nz, grid_type, 
        xvec, zvec, bin_rho, bin_vp, bin_vs);

    delete [] bin_rho;
    delete [] bin_vp;
    delete [] bin_vs;
    return 0;
}

int parameterization_bin_el_iso_loc(
    float *rho2d,
    float *lam2d,
    float *mu2d, 
    const float *x2d,
    const float *z2d,
    size_t nx, size_t nz,
    int grid_type,
    std::vector<float> &xvec, 
    std::vector<float> &zvec,
    float *bin_rho,
    float *bin_vp,
    float *bin_vs ) 
{

  size_t siz_line  = nx;
  size_t siz_slice = nx * nz;

  float slow_k = 1.0/(nz-1); // for print progress
  std::cout << "- discrete model from the binary file\n\n";
  for (size_t k = 0; k < nz; ++k) {
    printProgress(slow_k*k);
    for (size_t i = 0; i < nx; ++i) {
      size_t indx =  i + k * siz_line;
      size_t indx_z = indx;          // for vmap and curv: z
      size_t indx_x = i; // for vmap and cart: x, y 
      if (grid_type == GRID_CART) {
        indx_z = k;                // for cart: z
      } else if (grid_type == GRID_CURV) {
        indx_x = indx; // for curv: x, y
      }
      float rho, vp, vs;
      rho = TrilinearInterpolation(xvec, zvec, bin_rho, x2d[indx_x], z2d[indx_z]);    
      vp  = TrilinearInterpolation(xvec, zvec, bin_vp , x2d[indx_x], z2d[indx_z]);    
      vs  = TrilinearInterpolation(xvec, zvec, bin_vs , x2d[indx_x], z2d[indx_z]); 

      mu2d[indx]  = vs*vs*rho;
      lam2d[indx] = vp*vp*rho-2.0*mu2d[indx];
      rho2d[indx] = rho;
    }
  }

  return 0;
}
