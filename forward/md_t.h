#ifndef MD_EL_ISO_H
#define MD_EL_ISO_H

#include "gd_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int n1, n2, n3, n4;
  int nx, nz, ncmp;
  float *v3d; // allocated var

  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

  // flag to determine medium type
  int medium_type;
  int visco_type;

  // rho for all media
  float *rho;

  // for acustic
  float *kappa; // pointer to var

  // for isotropic media
  float *lambda; // pointer to var
  float *mu;

  // for visco attenuation
  float *Qs;

  // for anisotropic media
  float *c11;
  float *c13;
  float *c15;
  float *c33;
  float *c35;
  float *c55;

  float visco_Qs_freq;

} md_t;

/*************************************************
 * function prototype
 *************************************************/

int
md_init(gd_t *gd, md_t *md, int media_type, int visco_type);

int
md_import(md_t *md, char *in_dir);

int
md_export(gd_t *gd,
          md_t *md,
          char *output_dir);

int
md_gen_test_ac_iso(md_t *md);

int
md_gen_test_el_iso(md_t *md);

int
md_gen_test_el_vti(md_t *md);

int
md_gen_test_Qs(md_t *md, float Qs_freq);

int
md_gen_test_el_aniso(md_t *md);

int
md_rho_to_slow(float * rho, size_t siz_icmp);

#endif
