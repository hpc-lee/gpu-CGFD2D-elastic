#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include "alloc.h"
#include "cuda_common.h"


int init_gd_device(gd_t *gd, gd_t *gd_d)
{
  memcpy(gd_d,gd,sizeof(gd_t));
  return 0;
}

int init_md_device(md_t *md, md_t *md_d)
{
  size_t siz_icmp = md->siz_icmp;

  memcpy(md_d,md,sizeof(md_t));
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO)
  {
    md_d->rho    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->kappa  = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    CUDACHECK(cudaMemcpy(md_d->rho,    md->rho,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->kappa,  md->kappa,  sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
  }
  if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
  {
    md_d->rho    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->lambda = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->mu     = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    CUDACHECK(cudaMemcpy(md_d->rho,    md->rho,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->lambda, md->lambda, sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->mu,     md->mu,     sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
  }
  if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI)
  {
    md_d->rho    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c11    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c13    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c33    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c55    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    CUDACHECK(cudaMemcpy(md_d->rho,    md->rho,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c11,    md->c11,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c13,    md->c13,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c33,    md->c33,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c55,    md->c55,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
  }
  if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO)
  {
    md_d->rho    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c11    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c13    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c15    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c33    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c35    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    md_d->c55    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
    CUDACHECK(cudaMemcpy(md_d->rho,    md->rho,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c11,    md->c11,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c13,    md->c13,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c15,    md->c15,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c33,    md->c33,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c35,    md->c35,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(md_d->c55,    md->c55,    sizeof(float)*siz_icmp, cudaMemcpyHostToDevice));
  }

  return 0;
}

int init_fd_device(fd_t *fd, fd_wav_t *fd_wav_d)
{
  int max_len = fd->fdz_max_len; //=5 
  int max_lay = fd->num_of_fdz_op;
  fd_wav_d->fdz_len_d     = (int *) cuda_malloc(sizeof(int)*max_lay);
  fd_wav_d->fdx_coef_d    = (float *) cuda_malloc(sizeof(float)*max_len);
  fd_wav_d->fdz_coef_d    = (float *) cuda_malloc(sizeof(float)*max_len);
  fd_wav_d->fdz_coef_all_d    = (float *) cuda_malloc(sizeof(float)*max_len*max_lay);

  fd_wav_d->fdx_indx_d    = (int *) cuda_malloc(sizeof(int)*max_len);
  fd_wav_d->fdz_indx_d    = (int *) cuda_malloc(sizeof(int)*max_len);
  fd_wav_d->fdz_indx_all_d    = (int *) cuda_malloc(sizeof(int)*max_len*max_lay);

  fd_wav_d->fdx_shift_d    = (size_t *) cuda_malloc(sizeof(size_t)*max_len);
  fd_wav_d->fdz_shift_d    = (size_t *) cuda_malloc(sizeof(size_t)*max_len);
  fd_wav_d->fdz_shift_all_d    = (size_t *) cuda_malloc(sizeof(size_t)*max_len*max_lay);
  return 0;
}

int init_metric_device(gdcurv_metric_t *metric, gdcurv_metric_t *metric_d)
{
  size_t siz_icmp = metric->siz_icmp;

  memcpy(metric_d,metric,sizeof(gdcurv_metric_t));
  metric_d->jac     = (float *) cuda_malloc(sizeof(float)*siz_icmp);
  metric_d->xi_x    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
  metric_d->xi_z    = (float *) cuda_malloc(sizeof(float)*siz_icmp);
  metric_d->zeta_x   = (float *) cuda_malloc(sizeof(float)*siz_icmp);
  metric_d->zeta_z   = (float *) cuda_malloc(sizeof(float)*siz_icmp);

  CUDACHECK( cudaMemcpy(metric_d->jac,   metric->jac,   sizeof(float)*siz_icmp, cudaMemcpyHostToDevice) );
  CUDACHECK( cudaMemcpy(metric_d->xi_x,  metric->xi_x,  sizeof(float)*siz_icmp, cudaMemcpyHostToDevice) );
  CUDACHECK( cudaMemcpy(metric_d->xi_z,  metric->xi_z,  sizeof(float)*siz_icmp, cudaMemcpyHostToDevice) );
  CUDACHECK( cudaMemcpy(metric_d->zeta_x, metric->zeta_x, sizeof(float)*siz_icmp, cudaMemcpyHostToDevice) );
  CUDACHECK( cudaMemcpy(metric_d->zeta_z, metric->zeta_z, sizeof(float)*siz_icmp, cudaMemcpyHostToDevice) );
  return 0;
}

int init_src_device(src_t *src, src_t *src_d)
{
  int total_number = src->total_number;
  int max_ext      = src->max_ext;
  size_t temp_all     = (src->total_number) * (src->max_nt) * (src->max_stage);

  memcpy(src_d,src,sizeof(src_t));
  if(src->force_actived == 1) {
    src_d->Fx  = (float *) cuda_malloc(sizeof(float)*temp_all);
    src_d->Fz  = (float *) cuda_malloc(sizeof(float)*temp_all);
    CUDACHECK( cudaMemcpy(src_d->Fx,  src->Fx,  sizeof(float)*temp_all, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->Fz,  src->Fz,  sizeof(float)*temp_all, cudaMemcpyHostToDevice));
  }
  if(src->moment_actived == 1) {
    src_d->Mxx = (float *) cuda_malloc(sizeof(float)*temp_all);
    src_d->Mzz = (float *) cuda_malloc(sizeof(float)*temp_all);
    src_d->Mxz = (float *) cuda_malloc(sizeof(float)*temp_all);
    CUDACHECK( cudaMemcpy(src_d->Mxx, src->Mxx, sizeof(float)*temp_all, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->Mzz, src->Mzz, sizeof(float)*temp_all, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->Mxz, src->Mxz, sizeof(float)*temp_all, cudaMemcpyHostToDevice));
  }
  if(total_number>0)
  {
    src_d->ext_num  = (int *) cuda_malloc(sizeof(int)*total_number);
    src_d->it_begin = (int *) cuda_malloc(sizeof(int)*total_number);
    src_d->it_end   = (int *) cuda_malloc(sizeof(int)*total_number);
    src_d->ext_indx = (int *) cuda_malloc(sizeof(int)*total_number*max_ext);
    src_d->ext_coef = (float *) cuda_malloc(sizeof(float)*total_number*max_ext);

    CUDACHECK( cudaMemcpy(src_d->ext_num,  src->ext_num,  sizeof(int)*total_number, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->it_begin, src->it_begin, sizeof(int)*total_number, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->it_end,   src->it_end,   sizeof(int)*total_number, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->ext_indx, src->ext_indx, sizeof(int)*total_number*max_ext, cudaMemcpyHostToDevice));
    CUDACHECK( cudaMemcpy(src_d->ext_coef, src->ext_coef, sizeof(float)*total_number*max_ext, cudaMemcpyHostToDevice));
  }
  return 0;
}

int init_bdry_device(gd_t *gd, bdry_t *bdry, bdry_t *bdry_d)
{
  int nx = gd->nx;
  int nz = gd->nz;

  memcpy(bdry_d,bdry,sizeof(bdry_t));
  // copy bdryfree
  if (bdry_d->is_sides_free[CONST_NDIM-1][1] == 1)
  {
    bdry_d->vecVx2Vz2   = (float *) cuda_malloc(sizeof(float)*nx*CONST_NDIM*CONST_NDIM);

    CUDACHECK(cudaMemcpy(bdry_d->vecVx2Vz2, bdry->vecVx2Vz2, sizeof(float)*nx*CONST_NDIM*CONST_NDIM, cudaMemcpyHostToDevice));
  }

  // copy bdrypml
  if (bdry_d->is_enable_pml == 1)
  {
    for(int idim=0; idim<CONST_NDIM; idim++){
      for(int iside=0; iside<2; iside++){
        if(bdry_d->is_sides_pml[idim][iside] == 1){
          int npoints = bdry_d->num_of_layers[idim][iside] + 1;
          bdry_d->A[idim][iside]   = (float *) cuda_malloc(npoints * sizeof(float));
          bdry_d->B[idim][iside]   = (float *) cuda_malloc(npoints * sizeof(float));
          bdry_d->D[idim][iside]   = (float *) cuda_malloc(npoints * sizeof(float));
          CUDACHECK(cudaMemcpy(bdry_d->A[idim][iside],bdry->A[idim][iside],npoints*sizeof(float),cudaMemcpyHostToDevice));
          CUDACHECK(cudaMemcpy(bdry_d->B[idim][iside],bdry->B[idim][iside],npoints*sizeof(float),cudaMemcpyHostToDevice));
          CUDACHECK(cudaMemcpy(bdry_d->D[idim][iside],bdry->D[idim][iside],npoints*sizeof(float),cudaMemcpyHostToDevice));
          } else {
          bdry_d->A[idim][iside] = NULL;
          bdry_d->B[idim][iside] = NULL;
          bdry_d->D[idim][iside] = NULL;
        }
      }
    }

    for(int idim=0; idim<CONST_NDIM; idim++){
      for(int iside=0; iside<2; iside++){
        bdrypml_auxvar_t *auxvar_d = &(bdry_d->auxvar[idim][iside]);
        if(auxvar_d->siz_icmp > 0){
          auxvar_d->var = (float *) cuda_malloc(sizeof(float)*auxvar_d->siz_ilevel*auxvar_d->nlevel); 
          CUDACHECK(cudaMemset(auxvar_d->var,0,sizeof(float)*auxvar_d->siz_ilevel*auxvar_d->nlevel));
        } else {
        auxvar_d->var = NULL;
        }
      }
    }
  }
  // copy bdryexp
  if (bdry_d->is_enable_ablexp == 1)
  {
    bdry_d->ablexp_Ex = (float *) cuda_malloc(nx * sizeof(float));
    bdry_d->ablexp_Ez = (float *) cuda_malloc(nz * sizeof(float));
    CUDACHECK(cudaMemcpy(bdry_d->ablexp_Ex,bdry->ablexp_Ex,nx*sizeof(float),cudaMemcpyHostToDevice));
    CUDACHECK(cudaMemcpy(bdry_d->ablexp_Ez,bdry->ablexp_Ez,nz*sizeof(float),cudaMemcpyHostToDevice));
  }

  return 0;
}


int init_wave_device(wav_t *wav, wav_t *wav_d)
{
  size_t siz_ilevel = wav->siz_ilevel;
  int nlevel = wav->nlevel;
  memcpy(wav_d,wav,sizeof(wav_t));
  wav_d->v4d   = (float *) cuda_malloc(sizeof(float)*siz_ilevel*nlevel);
  CUDACHECK(cudaMemset(wav_d->v4d,0,sizeof(float)*siz_ilevel*nlevel));

  return 0;
}


int dealloc_gdcurv_device(gd_t gdcurv_d)
{
  CUDACHECK(cudaFree(gdcurv_d.x2d)); 
  CUDACHECK(cudaFree(gdcurv_d.z2d)); 

  return 0;
}

int dealloc_md_device(md_t md_d)
{
  if (md_d.medium_type == CONST_MEDIUM_ACOUSTIC_ISO)
  {
    CUDACHECK(cudaFree(md_d.rho   )); 
    CUDACHECK(cudaFree(md_d.kappa)); 
  }
  if (md_d.medium_type == CONST_MEDIUM_ELASTIC_ISO)
  {
    CUDACHECK(cudaFree(md_d.rho   )); 
    CUDACHECK(cudaFree(md_d.lambda)); 
    CUDACHECK(cudaFree(md_d.mu    )); 
  }
  if (md_d.medium_type == CONST_MEDIUM_ELASTIC_VTI)
  {
    CUDACHECK(cudaFree(md_d.rho)); 
    CUDACHECK(cudaFree(md_d.c11)); 
    CUDACHECK(cudaFree(md_d.c13)); 
    CUDACHECK(cudaFree(md_d.c33)); 
    CUDACHECK(cudaFree(md_d.c55)); 
  }
  if (md_d.medium_type == CONST_MEDIUM_ELASTIC_ANISO)
  {
    CUDACHECK(cudaFree(md_d.rho)); 
    CUDACHECK(cudaFree(md_d.c11)); 
    CUDACHECK(cudaFree(md_d.c13)); 
    CUDACHECK(cudaFree(md_d.c15)); 
    CUDACHECK(cudaFree(md_d.c33)); 
    CUDACHECK(cudaFree(md_d.c35)); 
    CUDACHECK(cudaFree(md_d.c55)); 
  }

  return 0;
}

int dealloc_fd_device(fd_wav_t fd_wav_d)
{
  CUDACHECK(cudaFree(fd_wav_d.fdz_len_d));

  CUDACHECK(cudaFree(fd_wav_d.fdx_coef_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_coef_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_coef_all_d));

  CUDACHECK(cudaFree(fd_wav_d.fdx_indx_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_indx_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_indx_all_d));

  CUDACHECK(cudaFree(fd_wav_d.fdx_shift_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_shift_d));
  CUDACHECK(cudaFree(fd_wav_d.fdz_shift_all_d));

  return 0;
}
int dealloc_metric_device(gdcurv_metric_t metric_d)
{
  CUDACHECK(cudaFree(metric_d.jac   )); 
  CUDACHECK(cudaFree(metric_d.xi_x  )); 
  CUDACHECK(cudaFree(metric_d.xi_z  )); 
  CUDACHECK(cudaFree(metric_d.zeta_x)); 
  CUDACHECK(cudaFree(metric_d.zeta_z)); 
  return 0;
}

int dealloc_src_device(src_t src_d)
{
  if(src_d.force_actived == 1)
  {
    CUDACHECK(cudaFree(src_d.Fx)); 
    CUDACHECK(cudaFree(src_d.Fz)); 
  }
  if(src_d.moment_actived == 1)
  {
    CUDACHECK(cudaFree(src_d.Mxx)); 
    CUDACHECK(cudaFree(src_d.Mzz)); 
    CUDACHECK(cudaFree(src_d.Mxz)); 
  }
  if(src_d.total_number > 0)
  {
    CUDACHECK(cudaFree(src_d.ext_num )); 
    CUDACHECK(cudaFree(src_d.ext_indx)); 
    CUDACHECK(cudaFree(src_d.ext_coef)); 
    CUDACHECK(cudaFree(src_d.it_begin)); 
    CUDACHECK(cudaFree(src_d.it_end  )); 
  }
  return 0;
}

int dealloc_bdry_device(bdry_t bdry_d)
{
  if (bdry_d.is_sides_free[CONST_NDIM-1][1] == 1)
  {
    CUDACHECK(cudaFree(bdry_d.vecVx2Vz2)); 
  }
  if (bdry_d.is_enable_pml == 1)
  {
    for(int idim=0; idim<CONST_NDIM; idim++){
      for(int iside=0; iside<2; iside++){
        if(bdry_d.is_sides_pml[idim][iside] == 1){
          CUDACHECK(cudaFree(bdry_d.A[idim][iside])); 
          CUDACHECK(cudaFree(bdry_d.B[idim][iside])); 
          CUDACHECK(cudaFree(bdry_d.D[idim][iside])); 
        }
      }
    }  
    for(int idim=0; idim<CONST_NDIM; idim++){
      for(int iside=0; iside<2; iside++){
        bdrypml_auxvar_t *auxvar_d = &(bdry_d.auxvar[idim][iside]);
        if(auxvar_d->siz_icmp > 0){
          CUDACHECK(cudaFree(auxvar_d->var)); 
        }
      }
    }  
  }
  if (bdry_d.is_enable_ablexp == 1)
  {
    CUDACHECK(cudaFree(bdry_d.ablexp_Ex)); 
    CUDACHECK(cudaFree(bdry_d.ablexp_Ez));
  }
  return 0;
}

int dealloc_wave_device(wav_t wav_d)
{
  CUDACHECK(cudaFree(wav_d.v4d)); 
  return 0;
}

