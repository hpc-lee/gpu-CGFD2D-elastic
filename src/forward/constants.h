#ifndef CONSTANTS_H
#define CONSTANTS_H

// consts
#define CONST_NDIM     2
#define CONST_NDIM_2   4 // 2 * ndim
#define CONST_TIJ_SIZE 3 // 
#define CONST_MAX_STRLEN 1024

#ifndef M_PI
#define PI 3.14159265358979323846264338327950288419716939937510
#else
#define PI M_PI
#endif

// medium type
#define CONST_MEDIUM_ACOUSTIC_ISO  1
#define CONST_MEDIUM_ELASTIC_ISO   2
#define CONST_MEDIUM_ELASTIC_VTI   3
#define CONST_MEDIUM_ELASTIC_ANISO 4

// visco type
#define CONST_VISCO_GRAVES_QS 1
#define CONST_VISCO_GMB 2

#define handle_nc_err(err)                       \
{                                                \
  if (err != NC_NOERR) {                         \
     fprintf(stderr,"nc error: %s\n", nc_strerror(err)); \
     exit(-1);                                   \
  }                                              \
}

#endif
