#ifndef MEDIA_DISCRETE_MODEL_H
#define MEDIA_DISCRETE_MODEL_H

// for C code call
#define MEDIA_USE_CART 1
#define MEDIA_USE_VMAP 2
#define MEDIA_USE_CURV 3

/*------------ bin2model --------------*/
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

/*--------------------------- layer2model --------------------- */
//---- 0. one component
int media_layer2model_onecmp(float *var2d,
                             const float *x2d, 
                             const float *z2d, 
                             size_t nx,
                             size_t nz,
                             int grid_type, 
                             const char *in_var_file,
                             const char *average_method);

//---- 1. elastic isotropic
int media_layer2model_ac_iso(
        float *rho2d,
        float *kappa2d,
        const float *x2d, 
        const float *z2d, 
        size_t nx,
        size_t nz,
        int grid_type, 
        const char *in_lay_file,
        const char *equivalent_medium_method);

//----  2. elastic isotropic
int media_layer2model_el_iso(
        float *lam2d,
        float *mu2d,
        float *rho2d,
        const float *x2d, 
        const float *z2d, 
        size_t nx,
        size_t nz,
        int grid_type, 
        const char *in_lay_file,
        const char *equivalent_medium_method);

//--- 3. elastic vti
int media_layer2model_el_vti(
        float *rho,
        float *c11,
        float *c33,
        float *c55,
        float *c13,
        const float *x2d,
        const float *z2d,
        size_t nx,
        size_t nz,
        int grid_type, 
        const char *in_lay_file, 
        const char *equivalent_medium_method);

//--- 4. elastic anisotropic/TTI
int media_layer2model_el_aniso(
        float *rho,
        float *c11, float *c13, float *c15,
        float *c33, float *c35, 
        float *c55,
        const float *x2d,
        const float *z2d,
        size_t nx,
        size_t nz,
        int grid_type, 
        const char *in_lay_file,
        const char *equivalent_medium_method) ; 

/*-------------- grid2model -------------*/

// --- 0. one component
int media_grid2model_onecmp(
    float *var2d,
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 1. acoustic isotropic
int media_grid2model_ac_iso(
    float *rho2d,
    float *kappa2d, 
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 2. elastic isotropic
int media_grid2model_el_iso(
    float *rho2d,
    float *lam2d,
    float *mu2d, 
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 3. elastic vti
int media_grid2model_el_vti(
    float *rho,
    float *c11, 
    float *c33,
    float *c55,
    float *c13,
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 4. elastic anisotropic/TTI
int media_grid2model_el_aniso(
    float *rho,
    float *c11, 
    float *c13,
    float *c15, 
    float *c33,
    float *c35, 
    float *c55,
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

#endif
