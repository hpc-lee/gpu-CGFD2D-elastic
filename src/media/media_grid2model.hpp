#ifndef _MEDIA_GRID2MODEL_
#define _MEDIA_GRID2MODEL_

#include "media_geometry2d.hpp"
#include "media_read_file.hpp"

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

int AssignGridMediaPara2Point(
    size_t ix, size_t iz, 
    Point2 A, 
    inter_t &interfaces,
    int media_type,
    std::vector<float> &var,
    std::vector<int> &NGz, int NL);

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for grid2model
void CalPointValue_grid(int media_type, 
                   inter_t &interfaces,
                   size_t inter_line,  
                   std::vector<float> &xvec,
                   int NI, 
                   std::set<int> &layer_indx,
                   Point2 &A,
                   std::vector<float> &elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var);


//- 0. assign the parameter directly (use the local values): one component 
void parametrization_grid_onecmp_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *var2d);

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_grid_ac_iso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho2d,
    float *kappa);

//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_grid_el_iso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho2d, 
    float *lam2d,
    float *mu2d);

//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_grid_el_vti_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_grid_el_aniso_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//======================= averaging/equivalent medium method =========================

/* 
 * For half-grid point, mark the interface number.  
 */
int LayerNumberAtPoint(
    Point2 A, 
    int NL,
    std::vector<int> &NGz, 
    inter_t &interfaces);

void MarkLayerNumber(
    int grid_type, 
    float *Hx, float *Hz,
    size_t nx, size_t nz,
    int NL, std::vector<int> &NGz,
    int *MaterNum,
    inter_t &interfaces);

//- 0.1 assign the parameter by volume harmonic averaging
//- one component
int parametrization_grid_onecmp_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var2d);

//- 0.2 assign the parameter by volume arithmetic averaging
//- one component
int parametrization_grid_onecmp_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var2d);

//- 1.1 assign the parameter by volume arithmetic and harmonic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho2d,
    float *kappa);

//- 1.2 assign the parameter by volume arithmetic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho2d,
    float *kappa);

//- 2.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic isotropic
// Moczo et al., 2002 
int parametrization_grid_el_iso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho2d,
    float *lam2d,
    float *mu2d );

//- 2.1 assign the parameter by volume arithmetic averaging method
//- elastic isotropic
int parametrization_grid_el_iso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho2d,
    float *lam2d,
    float *mu2d );

//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c13,
    float *rho);

//- 4.1 assign the parameter by volume arithmetic averaging method
//- elastic tti/anisotropic
int parametrization_grid_el_aniso_har(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//- 4.2 assign the parameter by volume arithmetic averaging method
//- elastic tti
int parametrization_grid_el_aniso_ari(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55,
    float *rho);

//- 2.3 Assign the parameter by tti equivalent medium parametrization method
//- elasic isotropic
int parametrization_grid_el_iso_tti(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35, 
    float *c55 ); 

void getInterfaceIntersection_grid(
    const std::vector<float> &x, float *z,                          
    const Point2 &v1, const Point2 &v2,    
    std::set<Point2> &intersectionList); 

int parametrization_grid_el_vti_tti(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho,
    float *c11,
    float *c13,
    float *c15,
    float *c33,
    float *c35,
    float *c55);


#endif // MEIDA_GRID2MODEL
