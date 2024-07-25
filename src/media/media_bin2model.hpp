#ifndef _MEDIA_BIN2MODEL_
#define _MEDIA_BIN2MODEL_

#include "media_geometry2d.hpp"
#include "media_read_file.hpp"

int media_bin2model_el_iso(
    float *rho2d,
    float *lam2d,
    float *mu2d, 
    const float *x2d,
    const float *z2d,
    size_t nx, size_t nz,
    float xmin, float xmax,
    int grid_type,
    int *bin_order,    // eg, [1, 0,]=[z, x] 0:x, 1:z
    int *bin_size,     // [ndim1, ndim2],
    float  *bin_spacing,  // [dh1, dh2],
    float  *bin_origin,   // [h0_1, h0_2],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs);

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
    float *bin_vs ); 

#endif  // MEDIA_BIN2MODEL 
