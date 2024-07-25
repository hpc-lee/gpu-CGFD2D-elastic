#ifndef GD_CURV_H
#define GD_CURV_H

#include "constants.h"

/*************************************************
 * structure
 *************************************************/

typedef enum {

  GD_TYPE_CART = 1,
  GD_TYPE_VMAP = 2,
  GD_TYPE_CURV = 3

} gd_type_t;

//  grid coordinate for both cart, vmap and curv
//    to reduce duplicated functions
typedef struct {

  gd_type_t type;

  int ni, nk;
  int nx, nz;
  int ni1, ni2;
  int nk1, nk2;
  int gni1, gni2;
  int gnk1, gnk2;

  int npoint_ghosts;
  int fdx_nghosts;
  int fdz_nghosts;

  // curvilinear coord name,
  char **index_name;

  int ncmp;
  float *v3d; // allocated var

  //to avoid ref x2d at different funcs
  float *x2d; // pointer to var
  float *z2d;

  // for cart grid
  float *x1d;
  float *z1d;
  float dx;
  float dz;
  // x0/y0/z0 for grid gen
  float x0_glob;
  float z0_glob;

  // min/max of this thread including ghost points
  float xmin, xmax;
  float zmin, zmax;

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

} gd_t;

//  for metric
typedef struct {
  int nx, nz, ncmp;
  float *v3d; // allocated var

  float *jac; // pointer to var
  float *xi_x;
  float *xi_z;
  float *zeta_x;
  float *zeta_z;

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gdcurv_metric_t;


/*************************************************
 * function prototype
 *************************************************/

int 
gd_curv_init(gd_t *gdcurv);

int 
gd_curv_metric_init(gd_t        *gd,
                    gdcurv_metric_t *metric);
int
gd_curv_metric_cal(gd_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int * fd_indx, float * fd_coef);

int 
mirror_symmetry(gd_t *gdcurv,float *v3d, int ncmp);

int 
geometric_symmetry(gd_t *gdcurv,float *v3d, int ncmp);

int
gd_curv_gen_cart(gd_t *gdcurv,
                 float dx, float x0,
                 float dz, float z0);

int
gd_curv_coord_import(gd_t *gdcurv, char *import_dir);

int
gd_curv_coord_export(gd_t *gdcurv,
                     char *output_dir);

int
gd_cart_coord_export(gd_t *gdcart,
                     char *output_dir);

int
gd_curv_metric_import(gd_t *gd, 
                      gdcurv_metric_t *metric, 
                      char *import_dir);

int
gd_curv_metric_export(gd_t *gd,
                      gdcurv_metric_t *metric,
                      char *output_dir);

int
gd_curv_set_minmax(gd_t *gdcurv);

int 
gd_cart_init_set(gd_t *gdcart,
                 float dx, float x0_glob,
                 float dz, float z0_glob);

int
gd_cart_coord_to_local_indx(gd_t *gdcart,
                            float sx,
                            float sz,
                            int   *ou_si, int *ou_sk,
                            float *ou_sx_inc, float *ou_sz_inc);

int
gd_curv_coord_to_local_indx(gd_t *gd,
                            float sx, float sz,
                            int *si, int *sk,
                            float *sx_inc, float *sz_inc);

int
gd_curv_depth_to_axis(gd_t *gd,
                      float sx,
                      float *sz);

float
linear_interp_1d(float ix, float *x, float *z);

int
gd_curv_coord2index_sample(float sx, float sz, 
    float *points_x, // x coord of all points
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_k,
    int    nx_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sk_curv);

float
gd_coord_get_x(gd_t *gd, int i, int k);

float
gd_coord_get_z(gd_t *gd, int i, int k);

int
gd_indx_set(gd_t *const gd,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_z,
            const int fdx_nghosts,
            const int fdz_nghosts,
            const int verbose);

int
gd_lindx_is_inner(int i, int k, gd_t *gd);

int
gd_pindx_is_inner(int i_phy, int k_phy, gd_t *gd);

int
gd_pindx_is_inner_i(int i_phy, gd_t *gd);

int
gd_pindx_is_inner_k(int k_phy, gd_t *gd);

int
gd_print(gd_t *gd);

#endif
