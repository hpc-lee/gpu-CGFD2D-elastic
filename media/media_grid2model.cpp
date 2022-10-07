/******************************************************************************
 *
 * This function is used to discretize a given grid model to a calculation grid.
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    09/2021: Created by Luqian Jiang 
 *
 *******************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>  
#include <numeric>
#include <set>
#include "media_geometry2d.hpp"
#include "media_utility.hpp"
#include "media_grid2model.hpp"
#include "media_read_file.hpp"

std::map<int, std::string> md2str = create_md2str_map();
extern int edgeTable[16];

/*============================ for C call ====================================*/
//--- 0. one component
int media_grid2model_onecmp(
    float *var2d,
    const float *x2d,
    const float *z2d,
    int nx,
    int nz,
    float Xmin, float Xmax,
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz;
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, NL, NGz, &interfaces);

    if (interfaces.media_type != ONE_COMPONENT) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_onecmp() only supports one_component,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);        
        fflush(stderr);
        exit(1);
     }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_onecmp_loc(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, var2d);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_onecmp_har(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, var2d);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_onecmp_ari(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, var2d);
    } else {  //default
        fprintf(stderr,"Error: Wrong average method %s for one_component.\n", 
                equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }

    fprintf(stdout, " - Done\n"); 
    return 0;
}

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
    const char *equivalent_medium_method)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, NL, NGz, &interfaces);

    if (interfaces.media_type != ACOUSTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_ac_iso() only supports acoustic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);        
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_ac_iso_loc(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, kappa2d);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_ac_iso_har(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, kappa2d);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_ac_iso_ari(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, kappa2d);
    } else { //default
        fprintf(stderr,"Error: Wrong average method %s for acoustic_isotropic media.\n", 
            equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }
    fprintf(stdout, " - Done\n"); 

    return 0;
}

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
    const char *equivalent_medium_method)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, NL, NGz, &interfaces);

    if (interfaces.media_type != ELASTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_el_iso() only supports elastic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);         
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_el_iso_loc(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, lam2d, mu2d);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_el_iso_har(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, lam2d, mu2d);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_el_iso_ari(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, rho2d, lam2d, mu2d);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_grid2model_el_iso(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_grid2model_el_aniso().\n", 
                       equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }
    fprintf(stdout, " - Done\n"); 

    return 0;
}


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
    const char *equivalent_medium_method)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, NL, NGz, &interfaces);
    int md_type = interfaces.media_type;

    if (md_type != ELASTIC_VTI_PREM && md_type != ELASTIC_VTI_THOMSEN && md_type != ELASTIC_VTI_CIJ) {
        fprintf(stderr, "Error: media_type=%s is not supported in media_grid2model_el_vti(),\n"\
                        "       it only supports elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);    
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_el_vti_loc(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_el_vti_har(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_el_vti_ari(nx, nz, x2d, z2d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_grid2model_el_vti(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_grid2model_el_aniso().\n", 
                       equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }
    fprintf(stdout, " - Done\n"); 

    return 0;
}

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
    const char *equivalent_medium_method)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    size_t siz_slice = nx*nz;

    // Read grid media file
    read_grid_file(in_media_file, Xmin, Xmax, NL, NGz, &interfaces);

    int md_type = interfaces.media_type;

    // the function does not support one component and acoustic wave
    if (md_type == ONE_COMPONENT){
        fprintf(stderr, "Error: media_type=one_component is not supported in media_grid2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_media_file);        
        fflush(stderr);
        exit(1);
    } else if (md_type == ACOUSTIC_ISOTROPIC){
        fprintf(stderr, "Error: media_type=acoustic_isotropic is not supported in media_grid2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_media_file);         
        fflush(stderr);
        exit(1);
    }

    //- isotropic: loc, har, ari, tti
    if (md_type == ELASTIC_ISOTROPIC) {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_iso_loc(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, rho, c13, c55);
            for (size_t i = 0; i < siz_slice; ++i) {
                c11[i] = c13[i] + 2.0*c55[i]; 
                c33[i] = c11[i]; 
                c15[i] = 0.0; 
                c35[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_iso_har(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces,rho, c13, c55);
            for (size_t i = 0; i < siz_slice; ++i) {
                c11[i] = c13[i] + 2.0*c55[i]; 
                c33[i] = c11[i]; 
                c15[i] = 0.0; 
                c35[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_iso_ari(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, rho, c13, c55);
            for (size_t i = 0; i < siz_slice; ++i) {
                c11[i] = c13[i] + 2.0*c55[i]; 
                c33[i] = c11[i]; 
                c15[i] = 0.0; 
                c35[i] = 0.0;
            }
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
            parametrization_grid_el_iso_tti(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, rho, c11, c13, c15, c33, c35, c55);
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s!\n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    // vti: loc, har, ari    
    } else if (md_type == ELASTIC_VTI_PREM || md_type == ELASTIC_VTI_THOMSEN || md_type == ELASTIC_VTI_CIJ) {

        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_vti_loc(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c13, rho);
            for (size_t i = 0; i < siz_slice; i++) { 
                c15[i] = 0.0; c35[i] = 0.0; 
            }

        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_vti_har(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c13, rho);
            for (size_t i = 0; i < siz_slice; i++) {
                c15[i] = 0.0; c35[i] = 0.0; 
            }

        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_vti_ari(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c13, rho);
            for (size_t i = 0; i < siz_slice; i++) {
                c15[i] = 0.0; c35[i] = 0.0; 
            }

        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
            parametrization_grid_el_iso_tti(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, rho, c11, c13, c15, c33, c35, c55);
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    } else {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_aniso_loc(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c13, c15, c33, c35, c55, rho);
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_aniso_har(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c13, c15, c33, c35, c55, rho);
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_aniso_ari(nx, nz, x2d, z2d, grid_type, 
                NL, NGz, interfaces, c11, c13, c15, c33, c35, c55, rho);
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium for tti
        } else { //default
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    }

    fprintf(stdout, " - Done\n"); 

    return 0;
}
/*=======================================================================================================*/

int AssignGridMediaPara2Point(
    size_t ix, size_t iz, 
    Point2 A, 
    inter_t &interfaces,
    int media_type,
    std::vector<float> &var,
    std::vector<int> &NGz, int NL) 
{
    /* All the given grid mesh is the same */
    size_t  NI = interfaces.NI;
    size_t  NX = interfaces.NX;
    float MINX = interfaces.MINX;
    float   DX = interfaces.DX;
    float MAXX = MINX + (NX-1)*DX; 
    int isPointOnInter = 0; 

    std::set<int> layer_indx;
    int indx = 0;
    for (int i = 0; i < NL-1; i++) {
        indx += NGz[i];
        layer_indx.insert(indx);
    }

    // if out of the MEDIA GRID MESH area, exit !
    PrintIsPointOutOfInterfaceRange(A, ix, iz, MINX, MAXX);

    std::vector<float> XVEC(NX);
    for (size_t i = 0; i < NX; i++)
        XVEC[i] = MINX + DX*i;

    /* 
     * For each interface, interpolate the elevation of the position
     *   to get where the Point A is located in xoy plane.
     */
    std::vector<float> elevation(NI, -FLT_MAX);
    for (int ni = 0; ni < NI; ni++) {
        /* Get the elevation for the corresponding x location */
        elevation[ni] = LinearInterpolation(XVEC, interfaces.elevation + ni*NX, A.x);
    }

    /* find which material is used, the z-grid is given from top to bottom */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    if (layer_indx.count(mi) && isEqual(A.z, elevation[mi])) 
        isPointOnInter = 1;

    CalPointValue_grid(media_type, interfaces, NX, XVEC, NI, layer_indx, A, elevation, mi, var);

    return isPointOnInter;
}


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
                   std::vector<float> &var)
{
    float dis_r = 1.0/2.0;
    if ( fabs(elevation[mi+1]-elevation[mi]) > 1e-6 ) {
        dis_r = (A.z - elevation[mi])/(elevation[mi+1] - elevation[mi]);
    }
    switch(media_type)
    {
    case ONE_COMPONENT: /* 0. var */
        /* If grid_z > elevation of top_interface, it given by the medium of top non-zero thickness layer */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.var + mi*inter_line, A.x); 
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.var + mi*inter_line, A.x); 
        }else { 
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float var0 = LinearInterpolation(xvec, interfaces.var +     mi*inter_line, A.x);
            float var1 = LinearInterpolation(xvec, interfaces.var + (mi+1)*inter_line, A.x);
            var[0]  = var0  + (var1-var0) * dis_r;
        }
    break;

    case ELASTIC_ISOTROPIC: /*1. rho, vp, vs*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs  + mi*inter_line, A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs  + mi*inter_line, A.x);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float vp0  = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            float vs0  = LinearInterpolation(xvec, interfaces.vs  + mi*inter_line, A.x);
            float rho0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            mi += 1;
            float vp1  = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            float vs1  = LinearInterpolation(xvec, interfaces.vs  + mi*inter_line, A.x);
            float rho1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = vp0  + (vp1-vp0)*dis_r;
            var[2] = vs0  + (vs1-vs0)*dis_r;
            var[0] = rho0 + (rho1-rho0)*dis_r;
        }
    break;

    case ELASTIC_VTI_PREM: /*2. rho, vph, vpv, vsv, eta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vph + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vpv + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.vsv + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.eta + mi*inter_line, A.x);     
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vph + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vpv + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.vsv + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.eta + mi*inter_line, A.x);
        }
        else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float vph0 = LinearInterpolation(xvec, interfaces.vph + mi*inter_line, A.x);
            float vpv0 = LinearInterpolation(xvec, interfaces.vpv + mi*inter_line, A.x);
            float vsv0 = LinearInterpolation(xvec, interfaces.vsv + mi*inter_line, A.x);
            float rho0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float eta0 = LinearInterpolation(xvec, interfaces.eta + mi*inter_line, A.x);
            mi += 1;
            float vph1 = LinearInterpolation(xvec, interfaces.vph + mi*inter_line, A.x);
            float vpv1 = LinearInterpolation(xvec, interfaces.vpv + mi*inter_line, A.x);
            float vsv1 = LinearInterpolation(xvec, interfaces.vsv + mi*inter_line, A.x);
            float rho1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float eta1 = LinearInterpolation(xvec, interfaces.eta + mi*inter_line, A.x);
            var[0] = rho0 + (rho1-rho0)*dis_r;
            var[1] = vph0 + (vph1-vph0)*dis_r;
            var[2] = vpv0 + (vpv1-vpv0)*dis_r;
            var[3] = vsv0 + (vsv1-vsv0)*dis_r;
            var[4] = eta0 + (eta1-eta0)*dis_r;
        }   
    break;

    case ELASTIC_VTI_THOMSEN: /*3. rho, vp0, vs0, epsilon, delta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho    + mi*inter_line , A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp0    + mi*inter_line , A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs0    + mi*inter_line , A.x);
            var[3] = LinearInterpolation(xvec, interfaces.epsilon+ mi*inter_line , A.x);
            var[4] = LinearInterpolation(xvec, interfaces.delta  + mi*inter_line , A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho    + mi*inter_line , A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp0    + mi*inter_line , A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs0    + mi*inter_line , A.x);
            var[3] = LinearInterpolation(xvec, interfaces.epsilon+ mi*inter_line , A.x);
            var[4] = LinearInterpolation(xvec, interfaces.delta  + mi*inter_line , A.x);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho0   = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line, A.x);
            float vp00   = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line, A.x);
            float vs00   = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line, A.x);
            float epsil0 = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line, A.x);  // epsilon
            float delta0 = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line, A.x);
            mi += 1;
            float rho1   = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line, A.x);
            float vp01   = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line, A.x);
            float vs01   = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line, A.x);
            float epsil1 = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line, A.x);  // epsilon
            float delta1 = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line, A.x);
            var[0] = rho0   + ( rho1   - rho0  )*dis_r;
            var[1] = vp00   + ( vp01   - vp00  )*dis_r;
            var[2] = vs00   + ( vs01   - vs00  )*dis_r;
            var[3] = epsil0 + ( epsil1 - epsil0)*dis_r;
            var[4] = delta0 + ( delta1 - delta0)*dis_r;
        }   
    break;

    case ELASTIC_VTI_CIJ: /*4. rho c11 c33 c55 c13 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);   
        } else  if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);   
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho_0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float c11_0 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c33_0 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c55_0 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float c13_0 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            mi += 1;
            float rho_1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float c11_1 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c33_1 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c55_1 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float c13_1 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = c11_0 + (c11_1 - c11_0)*dis_r;
            var[2] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[3] = c55_0 + (c55_1 - c55_0)*dis_r;
            var[4] = c13_0 + (c13_1 - c13_0)*dis_r;
        }   
    break;

    case ELASTIC_TTI_THOMSEN: /*5. rho, vp0, vs0, epsilon, delta, dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line , A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line , A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line , A.x);
            var[3] = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line , A.x);
            var[4] = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line , A.x);
            var[5] = LinearInterpolation(xvec, interfaces.dip     + mi*inter_line , A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line , A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line , A.x);
            var[2] = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line , A.x);
            var[3] = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line , A.x);
            var[4] = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line , A.x);
            var[5] = LinearInterpolation(xvec, interfaces.dip     + mi*inter_line , A.x);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho_0     = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line, A.x);
            float vp0_0     = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line, A.x);
            float vs0_0     = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line, A.x);
            float dip_0     = LinearInterpolation(xvec, interfaces.dip     + mi*inter_line, A.x);
            float epsil_0   = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line, A.x);
            float delta_0   = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line, A.x);
            mi += 1;
            float rho_1     = LinearInterpolation(xvec, interfaces.rho     + mi*inter_line, A.x);
            float vp0_1     = LinearInterpolation(xvec, interfaces.vp0     + mi*inter_line, A.x);
            float vs0_1     = LinearInterpolation(xvec, interfaces.vs0     + mi*inter_line, A.x);
            float dip_1     = LinearInterpolation(xvec, interfaces.dip     + mi*inter_line, A.x);
            float epsil_1   = LinearInterpolation(xvec, interfaces.epsilon + mi*inter_line, A.x);
            float delta_1   = LinearInterpolation(xvec, interfaces.delta   + mi*inter_line, A.x);

            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = vp0_0 + (vp0_1 - vp0_0)*dis_r;
            var[2] = vs0_0 + (vs0_1 - vs0_0)*dis_r;
            var[3] = epsil_0 + (epsil_1 - epsil_0)*dis_r;
            var[4] = delta_0 + (delta_1 - delta_0)*dis_r;
            var[5] = dip_0 + (dip_1 - dip_0)*dis_r;    
        }   
    break;

    case ELASTIC_ANISO_CIJ: /* 7. rho c11 c13 c15 c33 c35 c55 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c15 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[5] = LinearInterpolation(xvec, interfaces.c35 + mi*inter_line, A.x);
            var[6] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c15 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[5] = LinearInterpolation(xvec, interfaces.c35 + mi*inter_line, A.x);
            var[6] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float c11_0 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c13_0 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            float c15_0 = LinearInterpolation(xvec, interfaces.c15 + mi*inter_line, A.x);
            float c33_0 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c35_0 = LinearInterpolation(xvec, interfaces.c35 + mi*inter_line, A.x);
            float c55_0 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float rho_0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            mi += 1;
            float c11_1 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c13_1 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            float c15_1 = LinearInterpolation(xvec, interfaces.c15 + mi*inter_line, A.x);
            float c33_1 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c35_1 = LinearInterpolation(xvec, interfaces.c35 + mi*inter_line, A.x);
            float c55_1 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float rho_1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = c11_0 + (c11_1 - c11_0)*dis_r; 
            var[2] = c13_0 + (c13_1 - c13_0)*dis_r;
            var[3] = c15_0 + (c15_1 - c15_0)*dis_r;
            var[4] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[5] = c35_0 + (c35_1 - c35_0)*dis_r;
            var[6] = c55_0 + (c55_1 - c55_0)*dis_r;
        }     
    break;

    case ELASTIC_TTI_BOND: /* 6. rho c11 c33 c55 c13 dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            var[5] = LinearInterpolation(xvec, interfaces.dip + mi*inter_line, A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            var[2] = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            var[3] = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            var[4] = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            var[5] = LinearInterpolation(xvec, interfaces.dip + mi*inter_line, A.x);
        }
        else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }            float c11_0 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c33_0 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c55_0 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float c13_0 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            float rho_0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float dip_0 = LinearInterpolation(xvec, interfaces.dip + mi*inter_line, A.x);
            mi += 1;
            float c11_1 = LinearInterpolation(xvec, interfaces.c11 + mi*inter_line, A.x);
            float c33_1 = LinearInterpolation(xvec, interfaces.c33 + mi*inter_line, A.x);
            float c55_1 = LinearInterpolation(xvec, interfaces.c55 + mi*inter_line, A.x);
            float c13_1 = LinearInterpolation(xvec, interfaces.c13 + mi*inter_line, A.x);
            float rho_1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            float dip_1 = LinearInterpolation(xvec, interfaces.dip + mi*inter_line, A.x);
            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = c11_0 + (c11_1 - c11_0)*dis_r;
            var[2] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[3] = c55_0 + (c55_1 - c55_0)*dis_r;
            var[4] = c13_0 + (c13_1 - c13_0)*dis_r;
            var[5] = dip_0 + (dip_1 - dip_0)*dis_r;
        }   
    break;

    case ACOUSTIC_ISOTROPIC: /* 7. rho vp */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
        } else if (mi == NI-1) {
            var[0] = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[1] = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float vp0  = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            float rho0 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            mi+=1;
            float vp1  = LinearInterpolation(xvec, interfaces.vp  + mi*inter_line, A.x);
            float rho1 = LinearInterpolation(xvec, interfaces.rho + mi*inter_line, A.x);
            var[0] = rho0 + (rho1 - rho0) * dis_r;
            var[1] = vp0  + (vp1 - vp0 ) * dis_r;
        }
    break;

    default: // for self-check
        fprintf(stderr,"Error: Unknow media! (for code check, please contact Luqian Jiang)\n");
        fflush(stderr);
        exit(1);

    }     
} 

//- 0. assign the parameter directly (use the local values): one component 
void parametrization_grid_onecmp_loc(
    size_t nx, 
    size_t nz,
    const float *Gridx, 
    const float *Gridz,
    int grid_type,
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *var2d)
{
//    float slow_k = 1.0/(nz-1); // for print progress
//    std::cout << "- discrete model by local values:\n\n";
    for (size_t k = 0; k < nz; k++) {
        for (size_t i = 0; i < nx; i++) {
            std::vector<float> var(1, 0.0);

            size_t indx =  i + k * nx;
            size_t indx_z = indx;        // for vmap and curv: z
            size_t indx_x = i;           // for vmap and cart: x, y

            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }
            
            AssignGridMediaPara2Point(i, k,
                    Point2(Gridx[indx_x], Gridz[indx_z]), 
                    interfaces, ONE_COMPONENT, var, NGz, NL);
    
            var2d[indx] = var[0];     
        }
    }

}

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
    float *kappa)
{

//    float slow_k = 1.0/(nz-1); // for print progress
//    std::cout << "- discrete model by local values:\n\n";
    for (size_t k = 0; k < nz; ++k) {
        for (size_t i = 0; i < nx; ++i) {

            std::vector<float> var(2, 0.0); // rho, vp

            size_t indx =  i + k * nx;

            size_t indx_z = indx;          // for vmap and curv: z
            size_t indx_x = i;             // for vmap and cart: x, y

            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }

            AssignGridMediaPara2Point(i, k,
                    Point2(Gridx[indx_x], Gridz[indx_z]), 
                    interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);

            /*calculate output kappa */
            float vp = var[1], rho = var[0];
            kappa[indx] = vp*vp*rho;
            rho2d[indx] = rho; 
        }
    }
}


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
    float *mu2d)
{
    for (size_t k = 0; k < nz; ++k) {
        for (size_t i = 0; i < nx; ++i) {

            std::vector<float> var(3, 0.0); // rho, vp, vs

            size_t indx =  i + k*nx;

            size_t indx_z = indx;          // for vmap and curv: z
            size_t indx_x = i;             // for vmap and cart: x, y
            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }

            AssignGridMediaPara2Point(i, k,
                    Point2(Gridx[indx_x], Gridz[indx_z]), 
                    interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);

            /*calculate output lambda and mu */
            float vp = var[1], vs = var[2], rho = var[0];
            mu2d[indx]  = vs*vs*rho;
            lam2d[indx] = vp*vp*rho - 2.0*mu2d[indx];
            rho2d[indx] = rho; 
        }
    }
}


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
    float *rho)
{
    int media_type = interfaces.media_type;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t i = 0; i < nx; ++i) {

            size_t indx =  i + k * nx;

            size_t indx_z = indx;          // for vmap and curv: z
            size_t indx_x = i; // for vmap and cart: x, y
            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }

            std::vector<float> var(6, 0.0); 
            AssignGridMediaPara2Point(i, k,
                    Point2(Gridx[indx_x], Gridz[indx_z]), 
                    interfaces, media_type, var, NGz, NL);
            
            para2vti(var, media_type, 
                c11[indx], c33[indx], c55[indx], c13[indx], rho[indx]);
        }
    }
}

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
    float *rho)
{
    int media_type = interfaces.media_type;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t i = 0; i < nx; ++i) {

            size_t indx =  i + k*nx;

            size_t indx_z = indx;          // for vmap and curv: z
            size_t indx_x = i;             // for vmap and cart: x

            if (grid_type == GRID_CART) {
                indx_z = k;                // for cart: z
            } else if (grid_type == GRID_CURV) {
                indx_x = indx;
            }

            std::vector<float> var(7, 0.0); 

            AssignGridMediaPara2Point(i, k,
                        Point2(Gridx[indx_x], Gridz[indx_z]), 
                        interfaces, media_type, var, NGz, NL);

            para2tti(var, media_type, // return cij
                c11[indx], c13[indx], c15[indx],
                c33[indx], c35[indx], c55[indx], rho[indx]);
        }
        
    }
}


//======================= averaging/equivalent medium method =========================

/* 
 * For half-grid point, mark the interface number.  
 */
int LayerNumberAtPoint(
    Point2 A, 
    int NL,
    std::vector<int> &NGz, 
    inter_t &interfaces) 
{
    std::vector<float> elevation(NL, -FLT_MAX);
    float MINX = interfaces.MINX;
    float   DX = interfaces.DX;
    int     NX = interfaces.NX;
    float MAXX = MINX + (NX-1)*DX;  

    std::vector<float> XVEC(NX);
    for (int i = 0; i < NX; i++)
        XVEC[i] = MINX + DX*i;

    /* 
     * Mark the interfaces[0] and interfaces[accumulation of NGz[nl]], 
     * nl is form 0 to NL-2.
     */
    int ni = 0;
    for (int nl = 0; nl < NL; nl++) {
        /* Get the elevation for the corresponding location */
        elevation[nl] = LinearInterpolation(XVEC, interfaces.elevation + ni*NX, A.x);
        ni += NGz[nl];
    }

    /* Mark the layer number */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    /* If the elevation is above the surface, it is assigned by the first layer para. */
    if (mi == -1) 
        mi = findLastGreaterEqualIndex(elevation[0], elevation);

    return mi;
}

void MarkLayerNumber(
    int grid_type, 
    float *Hx, float *Hz,
    size_t nx, size_t nz,
    int NL, std::vector<int> &NGz,
    int *MaterNum,
    inter_t &interfaces) 
{
    /* 
     * Mark the layer number at the half grid,
     *  for the grid media format, just 0 - NL-1 are marked.
     */
    int siz_slice = nx*nz;
    if (grid_type == GRID_CART) {
        for (size_t k = 0; k < nz; k++) {
            for (size_t i = 0; i < nx; i++) {
                size_t indx =  i + k * nx; 
                MaterNum[indx] = LayerNumberAtPoint(Point2(Hx[i], Hz[k]), 
                                                    NL, NGz, interfaces);
            }
        }
    } else if (grid_type == GRID_VMAP) {
        for (size_t k = 0; k < nz; k++) {
            for (size_t i = 0; i < nx; i++) {
                size_t indx =  i + k * nx; 
                MaterNum[indx] = LayerNumberAtPoint(Point2(Hx[i], Hz[indx]), 
                                                    NL, NGz, interfaces);
            }
        }
    } else if (grid_type == GRID_CURV) {
        for (size_t i = 0; i < siz_slice; i++) {
            MaterNum[i] = LayerNumberAtPoint(Point2(Hx[i], Hz[i]), NL, NGz, interfaces);
        }
    }

}

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
    float *var2d) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_onecmp_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, var2d);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

             // There is more than one medium value in the half-grid mesh, subdivide the mesh.
            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float har_var = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {
                    std::vector<float> var(1, 0);
                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, ONE_COMPONENT, var, NGz, NL);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        har_var += (1.0/var[0]);
                    }
                }

                if (num_dis < 1) {
                    fprintf(stderr, "Too few sub-grid, please increase NSG!\n");
                    exit(1);
                }

                var2d[indx] = num_dis*1.0/har_var;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *var2d) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_onecmp_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, var2d);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i  + k * nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            // There is more than one medium value in the half-grid mesh, subdivide the mesh.
            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice,Hx, Hz);

                // recalculate the material value of the point
                float ari_var = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(1, 0);

                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, ONE_COMPONENT, var, NGz, NL);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        ari_var += (var[0]);
                    }
                }

                var2d[indx] = ari_var/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }
    
    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *kappa) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_ac_iso_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, rho2d, kappa);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
//    float slow_k = 1.0/(ny-2);
//    std::cout << "- equivalent medium parametrization:\n\n";

    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx = i + k*nx;

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                /* recalculate the material value of the point */
                float ari_rho   = 0.0;
                float har_kappa = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(3, 0);

                    int isPointOnInter = AssignGridMediaPara2Point(i, k, SubGrid[isg],
                        interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        float rho = var[0], vp = var[1];
                        ari_rho   += rho;
                        har_kappa += (1.0/(vp*vp*rho));
                    }
                }

                rho2d[indx] = ari_rho/num_dis;
                kappa[indx] = 1.0*num_dis/har_kappa;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }
    

    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *kappa) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_ac_iso_loc(nx,nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, rho2d, kappa);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                /* recalculate the material value of the point */
                float ari_rho   = 0.0;
                float ari_kappa = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(3, 0);

                    int isPointOnInter = AssignGridMediaPara2Point(i, k, SubGrid[isg],
                        interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);
                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        float rho = var[0], vp = var[1];
                        ari_rho   += rho;
                        ari_kappa += (vp*vp*rho);
                    }
                }

                rho2d[indx] = ari_rho/num_dis;
                kappa[indx] = ari_kappa/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}



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
    float *mu2d ) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_el_iso_loc(nx,nz, Gridx,Gridz,
        grid_type, NL, NGz, interfaces, rho2d, lam2d, mu2d);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {
                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                        nx, siz_slice, Hx, Hz);

                /* recalculate the material value of the point */
                float ari_rho   = 0.0;
                float har_kappa = 0.0;
                float har_mu    = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(3, 0);
                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);
                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        float rho = var[0], vp = var[1], vs = var[2];
                        float mu =  vs*vs*rho;
                        float lambda = vp*vp*rho - 2.0*mu;

                        ari_rho   += rho;
                        har_kappa += (1.0/(lambda + mu));
                        har_mu    += (1.0/mu);
                    }
                }

                har_mu    = num_dis*1.0/har_mu;
                har_kappa = num_dis*1.0/har_kappa;
                ari_rho   = ari_rho/num_dis; 

                lam2d[indx] = har_kappa - har_mu;
                mu2d[indx]  = har_mu;
                rho2d[indx] = ari_rho;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }
    


    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *mu2d ) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_el_iso_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, rho2d, lam2d, mu2d);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, 
            NL, NGz, MaterNum, interfaces);
    
    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float ari_rho = 0.0;
                float ari_lam = 0.0;
                float ari_mu  = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);

                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(3, 0.0);

                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        float vp = var[1], vs = var[2], rho = var[0];
                        float mu = vs*vs*rho;
                        float lambda = vp*vp*rho - 2.0*mu;

                        ari_rho += rho;
                        ari_lam += lambda;
                        ari_mu  += mu;
                    }
                }

                lam2d[indx] = ari_lam/(num_dis*1.0);
                mu2d[indx]  = ari_mu /(num_dis*1.0);
                rho2d[indx] = ari_rho/(num_dis*1.0);

                if(SubGrid != nullptr) delete []SubGrid;

            }

        }
    }
    
    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


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
    float *rho)
{
    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_vti_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);
   
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                        nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float ari_rho = 0.0;
                float har_c11 = 0.0;
                float har_c33 = 0.0;
                float har_c55 = 0.0;
                float har_c13 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(6, 0);

                    int isPointOnInter = AssignGridMediaPara2Point(i, k, SubGrid[isg],
                        interfaces, media_type, var, NGz, NL);
                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c33_p = 0.0;
                        float c55_p = 0.0;
                        float c13_p = 0.0, rho_p = 0.0;
    
                        para2vti(var, media_type, 
                            c11_p, c33_p, c55_p, c13_p, rho_p);
    
                        har_c11 += (1.0/c11_p);
                        har_c13 += (1.0/c13_p);
                        har_c33 += (1.0/c33_p);
                        har_c55 += (1.0/c55_p);
                        ari_rho += rho_p;
                    }   
                }

                c11[indx] = 1.0*num_dis/har_c11;
                c13[indx] = 1.0*num_dis/har_c13;
                c33[indx] = 1.0*num_dis/har_c33;
                c55[indx] = 1.0*num_dis/har_c55;
                rho[indx] = ari_rho/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


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
    float *rho)
{
    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_vti_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx,&Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            /* 
             * There is more than one medium value in the half-grid mesh, 
             *  subdivide the mesh.
             */
            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                // recalculate the material value of the point
                float ari_rho = 0.0;
                float ari_c11 = 0.0;
                float ari_c33 = 0.0;
                float ari_c55 = 0.0;
                float ari_c13 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(6, 0);

                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, media_type, var, NGz, NL);
                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c33_p = 0.0;
                        float c55_p = 0.0;
                        float c13_p = 0.0, rho_p = 0.0;
                        para2vti(var, media_type, 
                            c11_p, c33_p, c55_p, c13_p, rho_p);
    
                        ari_c11 += c11_p;
                        ari_c13 += c13_p;
                        ari_c33 += c33_p;
                        ari_c55 += c55_p;
                        ari_rho += rho_p;
                    }
                }

                c11[indx] = ari_c11/num_dis;
                c13[indx] = ari_c13/num_dis;
                c33[indx] = ari_c33/num_dis;
                c55[indx] = ari_c55/num_dis;
                rho[indx] = ari_rho/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *rho)
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_aniso_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, c11, c13, c15, c33, c35, c55, rho);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 
            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                    nx, siz_slice, Hx, Hz);

                /// recalculate the material value of the point
                float ari_rho = 0.0;
                float har_c11 = 0.0;
                float har_c13 = 0.0;
                float har_c15 = 0.0;
                float har_c33 = 0.0;
                float har_c35 = 0.0;
                float har_c55 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(8, 0.0);

                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, media_type, var, NGz, NL);

                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c13_p = 0.0;
                        float c15_p = 0.0, c33_p = 0.0;
                        float c35_p = 0.0, c55_p = 0.0;
                        float rho_p = 0.0;
                        para2tti(var, media_type, 
                            c11_p, c13_p, c15_p, 
                            c33_p, c35_p, c55_p, rho_p);
    
                        har_c11 += (1.0/c11_p);                
                        har_c13 += (1.0/c13_p);
                        har_c15 += (1.0/c15_p);
                        har_c33 += (1.0/c33_p);
                        har_c35 += (1.0/c35_p);
                        har_c55 += (1.0/c55_p);
                        ari_rho += rho_p;
                    }
                }

                c11[indx] = (1.0*num_dis)/har_c11;
                c13[indx] = (1.0*num_dis)/har_c13;
                c15[indx] = (1.0*num_dis)/har_c15;
                c33[indx] = (1.0*num_dis)/har_c33;
                c35[indx] = (1.0*num_dis)/har_c35;
                c55[indx] = (1.0*num_dis)/har_c55;
                
                rho[indx] = ari_rho/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }

        }
    }
    

    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

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
    float *rho)
{
    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_aniso_loc(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, c11, c13, c15, c33, c35, c55, rho);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k*nx; 

            /* Check if the corresponding the half-grid mesh have different values */
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            /* 
             * There is more than one medium value in the half-grid mesh, 
             *  subdivide the mesh.
             */
            if ( NumOfValues(v, NL) > 1) {

                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                            nx, siz_slice, Hx, Hz);

                /// recalculate the material value of the point
                float ari_rho = 0.0;
                float ari_c11 = 0.0;
                float ari_c13 = 0.0;
                float ari_c15 = 0.0;
                float ari_c33 = 0.0;
                float ari_c35 = 0.0;
                float ari_c55 = 0.0;

                Point2 *SubGrid = MeshSubdivide(M);
                int nsg = (NG+1)*(NG+1);
                int num_dis = nsg;
                for (int isg = 0; isg < nsg; isg++) {

                    std::vector<float> var(8, 0.0);

                    int isPointOnInter = AssignGridMediaPara2Point(
                        i, k, SubGrid[isg], interfaces, media_type, var, NGz, NL);
                    if (isPointOnInter == 1) {
                        num_dis--;
                    } else{
                        // for every sub-point, transfer para to cij
                        float c11_p = 0.0, c13_p = 0.0;
                        float c15_p = 0.0, c33_p = 0.0;
                        float c35_p = 0.0, c55_p = 0.0;
                        float rho_p = 0.0;
                        para2tti(var, media_type, 
                            c11_p, c13_p, c15_p, 
                            c33_p, c35_p, c55_p, rho_p);

                        ari_c11 += c11_p;
                        ari_c13 += c13_p;
                        ari_c15 += c15_p;
                        ari_c33 += c33_p;
                        ari_c35 += c35_p;
                        ari_c55 += c55_p;
                        ari_rho += rho_p;
                    }
                }

                c11[indx] = ari_c11/num_dis;
                c13[indx] = ari_c13/num_dis;
                c15[indx] = ari_c15/num_dis;
                c33[indx] = ari_c33/num_dis;
                c35[indx] = ari_c35/num_dis;
                c55[indx] = ari_c55/num_dis;
                
                rho[indx] = ari_rho/num_dis;
 
                if (SubGrid != nullptr) delete[] SubGrid;
            }
        } 
    }

    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


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
    float *c55 ) 
{

    size_t siz_slice = nz * nx;
    int NI = interfaces.NI;
    size_t  NX = interfaces.NX;
    float MINX = interfaces.MINX;
    float   DX = interfaces.DX;

    std::vector<float> xvec(NX);
    for (int i = 0; i < NX; i++)
        xvec[i] = MINX + DX*i;

    /* assign the local value first.*/
    parametrization_grid_el_iso_har(nx, nz, Gridx,Gridz,
        grid_type, NL, NGz, interfaces, rho, c13, c55);
    for (int i = 0; i < siz_slice; i++) {
        c11[i] = c13[i] + 2.0*c55[i];
        c33[i] = c11[i];
        c15[i] = 0.0; c35[i] = 0.0;
    }
    // if layer < 2, no need tti equivalent and return
    if (NL < 2) {
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);


    /* tti equivalent medium */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) == 2) {
                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                        nx, siz_slice, Hx, Hz);

                //- convert value to sign +(>= 0), -(< 0)
                int sum = std::accumulate(v.begin(), v.end(), 0);
                int interNum = -1;
                sum = sum/4+1; 
                for (auto &sign:v) {
                    interNum = std::max(interNum, sign);
                    sign -= sum;
                }

                int cubeState = 0;
                if (v[0] < 0) cubeState |= 1;
                if (v[1] < 0) cubeState |= 2;
                if (v[2] < 0) cubeState |= 4;
                if (v[3] < 0) cubeState |= 8;


                // get intersection points, use set to avoid duplication
                std::set<Point2> intersectionList;

                if (edgeTable[cubeState] & 1)  // v[0] v[1]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[0], M.v[1], intersectionList);

                if (edgeTable[cubeState] & 2) // v[1] v[2]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[1], M.v[2], intersectionList);

                if (edgeTable[cubeState] & 4) // v[2] v[3]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[2], M.v[3], intersectionList);

                if (edgeTable[cubeState] & 8) // v[3] v[0]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[3], M.v[0], intersectionList);                

                // tti equivalent
                if (intersectionList.size() == 2) {
                    std::set<Point2>::iterator it = intersectionList.begin();
                    Point2 P1 = *it;
                    it++;
                    Point2 P2 = *it;
                    Vector2 x_new = P2 - P1;
                    if (P1.x > P2.x) 
                        x_new = P1-P2;

                    // -pi - pi
                    double theta = atan2(x_new.z,x_new.x);

                    // for tti
                    float tmp1 = 0.0; // lam/lam2mu 
                    float tmp2 = 0.0; // lam2mu-lam^2/lam2mu
                    float harM = 0.0;
                    float harmu = 0.0;

                    Point2 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0.0);
                        int isPointOnInter = AssignGridMediaPara2Point(
                            i, k, SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float vp = var[1], vs = var[2], rho = var[0];
                            float mu = vs*vs*rho, lam2mu = vp*vp*rho;
                            float lam = lam2mu-2.0*mu;
                            tmp1 += (lam/lam2mu);
                            tmp2 += (lam2mu-lam*lam/lam2mu);
                            harM += (1/lam2mu);
                            harmu+= (1/mu);
                        }
                    }

                    tmp1 /= num_dis;
                    tmp2 /= num_dis;
                    harM  = num_dis*1.0/harM;
                    harmu = num_dis*1.0/harmu;

                    float A = tmp2 + harM*tmp1*tmp1;
                    float B = harM*tmp1;

                    BondTransform2d(A, B, 0, harM, 0, harmu, theta,
                        c11[indx], c13[indx], c15[indx], c33[indx], c35[indx], c55[indx]);
                }

            } // numofvalue == 2
        }
    }
    
    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


void getInterfaceIntersection_grid(
    const std::vector<float> &x, float *z,                              // for interface
    const Point2 &v1, const Point2 &v2,    // the edge
    std::set<Point2> &intersectionList)    // segment(vertex1, vertex2)
{
    size_t npoint = x.size();
    for (int i = 0; i < npoint-1; i++) {
        Point2 b1(x[i], z[i]);
        Point2 b2(x[i+1], z[i+1]);
        if (SegmentProperIntersection(v1, v2, b1, b2)) {
            Point2 intersectionP = GetLineIntersection(v1, v2-v1, b1, b2-b1);
            intersectionList.insert(intersectionP);
        } else{
            if (b1 == v1 || b1 == v2) {
                intersectionList.insert(b1);
            } else if (b2 == v1 || b2 == v2)  {
                intersectionList.insert(b2); 
            } else{
                if (isPointOnSegment(b1, v1, v2)) intersectionList.insert(b1);
                if (isPointOnSegment(b2, v1, v2)) intersectionList.insert(b2);
                if (isPointOnSegment(v1, b1, b2)) intersectionList.insert(v1);
                if (isPointOnSegment(v2, b1, b2)) intersectionList.insert(v2);
            }
        }
    }
}

//- 3.3 assign the parameter by tti arithmetic and harmonic averaging method
//- elastic vti
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
    float *c55)
{
    size_t siz_slice = nz * nx;
    int media_type = interfaces.media_type;
    size_t  NX = interfaces.NX;
    float MINX = interfaces.MINX;
    float   DX = interfaces.DX;

    std::vector<float> xvec(NX);
    for (int i = 0; i < NX; i++)
        xvec[i] = MINX + DX*i;

    /* assign the local value first.*/
    parametrization_grid_el_vti_ari(nx, nz, Gridx, Gridz,
        grid_type, NL, NGz, interfaces, c11, c33, c55, c13, rho);

    for (int i = 0; i < siz_slice; i++) {
        c15[i] = 0.0; c35[i] = 0.0;
    }

    if (NL < 2) {
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, nz, grid_type, Gridx, Gridz, &Hx, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_slice];
    MarkLayerNumber(grid_type, Hx, Hz, nx, nz, NL, NGz, MaterNum, interfaces);


    /* tti equivalent medium */ 
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t i = 1; i < nx-1; i++) {
            size_t indx =  i + k * nx; 
            std::vector<int> v(4);
            v[0] = MaterNum[indx-1-nx];
            v[1] = MaterNum[indx  -nx];
            v[2] = MaterNum[indx     ];
            v[3] = MaterNum[indx-1   ];

            if ( NumOfValues(v, NL) == 2) {
                Mesh2 M = GenerateHalfMesh(grid_type, i, k, indx, 
                        nx, siz_slice, Hx, Hz);

                //- convert value to sign +(>= 0), -(< 0)
                int sum = std::accumulate(v.begin(), v.end(), 0);
                int interNum = -1;
                sum = sum/4+1; 
                for (auto &sign:v) {
                    interNum = std::max(interNum, sign);
                    sign -= sum;
                }

                int cubeState = 0;
                if (v[0] < 0) cubeState |= 1;
                if (v[1] < 0) cubeState |= 2;
                if (v[2] < 0) cubeState |= 4;
                if (v[3] < 0) cubeState |= 8;


                // get intersection points, use set to avoid duplication
                std::set<Point2> intersectionList;

                if (edgeTable[cubeState] & 1)  // v[0] v[1]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[0], M.v[1], intersectionList);

                if (edgeTable[cubeState] & 2) // v[1] v[2]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[1], M.v[2], intersectionList);

                if (edgeTable[cubeState] & 4) // v[2] v[3]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[2], M.v[3], intersectionList);

                if (edgeTable[cubeState] & 8) // v[3] v[0]
                    getInterfaceIntersection_grid(xvec, 
                        interfaces.elevation + interNum*NX,
                        M.v[3], M.v[0], intersectionList);                

                // tti equivalent
                if (intersectionList.size() == 2) {
                    std::set<Point2>::iterator it = intersectionList.begin();
                    Point2 P1 = *it;
                    it++;
                    Point2 P2 = *it;
                    Vector2 x_new = P2 - P1;
                    if (P1.x > P2.x) 
                        x_new = P1-P2;

                    // -pi - pi
                    double theta = atan2(x_new.z,x_new.x);

                    // for tti
                    float tmp1 = 0.0; // lam/lam2mu 
                    float tmp2 = 0.0; // lam2mu-lam^2/lam2mu
                    float har33 = 0.0;
                    float har55 = 0.0;

                    Point2 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(6, 0.0);
                        int isPointOnInter = AssignGridMediaPara2Point(
                            i, k, SubGrid[isg], interfaces, media_type, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float c11_p = 0.0, c33_p = 0.0;
                            float c55_p = 0.0;
                            float c13_p = 0.0, rho_p = 0.0;
    
                            para2vti(var, media_type, 
                                c11_p, c33_p, c55_p, c13_p, rho_p);
    
                            tmp1 += (c13_p/c33_p);
                            tmp2 += (c11_p-c13_p*c13_p/c33_p);
                            har33 += (1/c33_p);
                            har55 += (1/c55_p);
                        }
                    }

                    tmp1 /= num_dis;
                    tmp2 /= num_dis;
                    har33 = num_dis*1.0/har33;
                    har55 = num_dis*1.0/har55;

                    float c11_z = tmp2 + har33*tmp1*tmp1;
                    float c13_z = har33*tmp1;
                    float c33_z = har33;
                    float c55_z = har55;

                    BondTransform2d(c11_z, c13_z, 0, c33_z, 0, c55_z, theta,
                        c11[indx], c13[indx], c15[indx], c33[indx], c35[indx], c55[indx]);
                }

            } // numofvalue == 2
        }
    }

    if (Hx != nullptr) delete[] Hx;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}
