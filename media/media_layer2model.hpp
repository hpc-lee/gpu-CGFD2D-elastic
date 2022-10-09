#ifndef _MEDIA_LAYER2MODEL_
#define _MEDIA_LAYER2MODEL_

#include "media_geometry2d.hpp"
#include "media_utility.hpp"

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
        const char *equivalent_medium_method); 


int AssignLayerMediaPara2Point(
    size_t ix, size_t iz,         /* To print Error messages */ 
    Point2 A,  
    inter_t *interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var); 

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for layer2model
void CalPointValue_layer(int media_type, 
                   inter_t *interfaces,
                   Point2 &A,
                   std::vector<float> &elevation,  /*the elevation of point A at the projection position of the interface mesh. */
                   std::vector<int> &internum4elev,
                   int mi,
                   std::vector<float> &var);

//- 0. assign the parameter directly (use the local values): one component 
void parametrization_layer_onecmp_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *var2d);

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_layer_ac_iso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type, 
    inter_t *interfaces,
    float *kappa,
    float *rho2d);

//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_layer_el_iso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *lam2d,
    float *mu2d,
    float *rho2d);

//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_layer_el_vti_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_layer_el_aniso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

void MarkInterfaceNumber(
        int grid_type,
        float *Hx, float *Hz,
        size_t nx, size_t nz,
        int *MaterNum, // nx*ny*nz
        inter_t *interfaces);

//- 0.1 Assign the parameter by volume harmonic averaging
//- one component
void parametrization_layer_onecmp_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *var2d);

//- 0.1 Assign the parameter by volume arithmetic averaging
//- one component
void parametrization_layer_onecmp_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *var2d);

//- 1.0 Assign the parameter by volume harmonic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type,
    inter_t *interfaces,
    float *kappa, 
    float *rho2d);

//- 1.1 Assign the parameter by volume arithmetic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type,
    inter_t *interfaces,
    float *kappa, 
    float *rho2d) ;

//- 2.1 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic isotropic
// ref: Moczo, 2002, 3D Heterogeneous Staggered-Grid Finite-Difference Modeling of Seismic
//      Motion with Volume Harmonic and Arithmetic Averaging of Elastic Moduli and Densities
void parametrization_layer_el_iso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *lam2d,
    float *mu2d,
    float *rho2d) ;

//- 2.2 Assign the parameter by volume arithmetic averaging method 
//- elasic isotropic
void parametrization_layer_el_iso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type,
    inter_t *interfaces,
    float *lam2d,
    float *mu2d,
    float *rho2d);

//- 3.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic vti
void parametrization_layer_el_vti_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type, 
    inter_t *interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 3.1 Assign the parameter by volume arithmetic averaging method 
//- elasic vti
void parametrization_layer_el_vti_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 4.1 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic aniso
void parametrization_layer_el_aniso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//- 4.2 Assign the parameter by volume arithmetic averaging method 
//- elasic tti
void parametrization_layer_el_aniso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int media_type,
    inter_t *interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);


//- 2.3 Assign the parameter by tti equivalent medium parametrization method
//- elasic isotropic
void parametrization_layer_el_iso_tti(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *rho,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55); 

void getInterfaceIntersection_layer(
    int npoint, float *x, float *z,   // for interface
    const Point2 &v1, const Point2 &v2,
    std::set<Point2> &intersectionList);

void parametrization_layer_el_vti_tti(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int media_type, 
    inter_t *interfaces,
    float *rho,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55);


#endif /* __MEDID_LAYER2MODEL__ */
