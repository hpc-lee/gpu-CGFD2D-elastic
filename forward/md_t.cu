/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "constants.h"
#include "fdlib_mem.h"
#include "md_t.h"

int
md_init(gd_t *gd, md_t *md, int media_type, int visco_type)
{
  int ierr = 0;

  md->nx   = gd->nx;
  md->nz   = gd->nz;

  md->siz_iz   = md->nx;
  md->siz_icmp  = md->nx * md->nz;

  // media type
  md->medium_type = media_type;

  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO)
  {
    md->ncmp = 2; // rho + kappa
  }
  else if (media_type == CONST_MEDIUM_ELASTIC_ISO)
  {
      md->ncmp = 3; // rho + lambda + mu
  } else if (media_type == CONST_MEDIUM_ELASTIC_VTI) {
      md->ncmp = 5; // c11 13 33 55 + rho
  } else if (media_type == CONST_MEDIUM_ELASTIC_ANISO){
      md->ncmp = 7; // 11, 13, 15, 33, 35, 55, rho
  } else{
      fprintf(stderr,"ERROR: media_type=%d is not implemented\n",media_type);
      exit(1);
  }

  // visco
  md->visco_type = visco_type;
  if (visco_type == CONST_VISCO_GRAVES_QS) {
   md->ncmp += 1;
  } 
 // else {
 //     fprintf(stderr,"ERROR: visco_type=%d is not implemented\n",visco_type);
 //     exit(1);
 // }

  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */
  
  // vars
  md->v3d = (float *) fdlib_mem_calloc_1d_float(
                          md->siz_icmp * md->ncmp,
                          0.0, "md_init");
  if (md->v3d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(md->ncmp,
                                                         0,
                                                         "medium_init");

  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(md->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "medium_init");

  // set pos
  for (int icmp=0; icmp < md->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * md->siz_icmp;
  }

  // init
  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","rho");
  md->rho = md->v3d + cmp_pos[icmp];

  // acoustic iso
  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","kappa");
    md->kappa = md->v3d + cmp_pos[icmp];
  }

  // iso
  if (media_type == CONST_MEDIUM_ELASTIC_ISO) 
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","lambda");
    md->lambda = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","mu");
    md->mu = md->v3d + cmp_pos[icmp];
  }

  // vti
  if (media_type == CONST_MEDIUM_ELASTIC_VTI)
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v3d + cmp_pos[icmp];
  }

  // aniso
  if (media_type == CONST_MEDIUM_ELASTIC_ANISO)
  {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c15");
    md->c15 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c35");
    md->c35 = md->v3d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v3d + cmp_pos[icmp];
  }

  // plus Qs
  if (visco_type == CONST_VISCO_GRAVES_QS) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","Qs");
    md->Qs = md->v3d + cmp_pos[icmp];
  }
  
  // set pointer
  md->cmp_pos  = cmp_pos;
  md->cmp_name = cmp_name;

  return ierr;
}

//
//
//

int
md_import(md_t *md, char *in_dir)
{
  int ierr = 0;

  char in_file[CONST_MAX_STRLEN];
  
  int ncid, varid;
  
  // construct file name
  sprintf(in_file, "%s/media.nc", in_dir);
  
  // read in nc
  ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  handle_nc_err(ierr);
  
  for (int icmp=0; icmp < md->ncmp; icmp++) {
      ierr = nc_inq_varid(ncid, md->cmp_name[icmp], &varid);
      handle_nc_err(ierr);
  
      ierr = nc_get_var_float(ncid,varid,md->v3d + md->cmp_pos[icmp]);
      handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);  handle_nc_err(ierr);

  return ierr;
}

int
md_export(gd_t  *gd,
          md_t  *md,
          char *output_dir)
{
  int ierr = 0;

  size_t * m3d_pos   = md->cmp_pos;
  char  ** m3d_name  = md->cmp_name;
  int  number_of_vars = md->ncmp;
  int  nx = md->nx;
  int  nz = md->nz;
  int  ni1 = gd->ni1;
  int  nk1 = gd->nk1;
  int  ni  = gd->ni;
  int  nk  = gd->nk;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/media.nc", output_dir);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  handle_nc_err(ierr);

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[1]);
  handle_nc_err(ierr);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);
  handle_nc_err(ierr);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, m3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
    handle_nc_err(ierr);
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
  handle_nc_err(ierr);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = md->v3d + m3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
    handle_nc_err(ierr);
  }
  
  // close file
  ierr = nc_close(ncid);
  handle_nc_err(ierr);

  return ierr;
}

/*
 * test
 */

int
md_gen_test_ac_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *kappa3d = md->kappa;
  float *rho3d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        float Vp=3000.0;
        float rho=1500.0;
        float kappa = Vp*Vp*rho;
        kappa3d[iptr] = kappa;
        rho3d[iptr] = rho;
      }
  }

  return ierr;
}

int
md_gen_test_el_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  float *lam3d = md->lambda;
  float  *mu3d = md->mu;
  float *rho3d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        float Vp=3000.0;
        float Vs=2000.0;
        float rho=1500.0;
        float mu = Vs*Vs*rho;
        float lam = Vp*Vp*rho - 2.0*mu;
        lam3d[iptr] = lam;
         mu3d[iptr] = mu;
        rho3d[iptr] = rho;
      }
  }

  return ierr;
}

int
md_gen_test_el_vti(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float rho=1500.0;

        md->rho[iptr] = rho;

	      md->c11[iptr] = 25.2*1e9;//lam + 2.0f*mu;
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c55[iptr] = 5.12*1e9;//mu;
        //-- Vp ~ sqrt(c11/rho) = 4098
      }
  }

  return ierr;
}

int
md_gen_test_el_aniso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;

        float rho=1500.0;

        md->rho[iptr] = rho;

	      md->c11[iptr] = 25.2*1e9;//lam + 2.0f*mu;
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c15[iptr] = 0.0;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c35[iptr] = 0.0;
	      md->c55[iptr] = 5.12*1e9;//mu;

        //-- Vp ~ sqrt(c11/rho) = 4098
      }
  }

  return ierr;
}

int
md_gen_test_Qs(md_t *md, float Qs_freq)
{
  int ierr = 0;

  int nx = md->nx;
  int nz = md->nz;
  int siz_iz = md->siz_iz;

  md->visco_Qs_freq = Qs_freq;

  float *Qs = md->Qs;

  for (size_t k=0; k<nz; k++)
  {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + k * siz_iz;
        Qs[iptr] = 20;
      }
  }

  return ierr;
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

int
md_rho_to_slow(float * rho, size_t siz_volume)
{
  int ierr = 0;

  for (size_t iptr=0; iptr<siz_volume; iptr++) {
    if (rho[iptr] > 1e-10) {
      rho[iptr] = 1.0 / rho[iptr];
    } else {
      rho[iptr] = 0.0;
    }
  }

  return ierr;
}
