/*******************************************************************************
 * Curvilinear Grid Finite Difference 2D Wave Propagation Simulation 
 *
 * Copyright (c) 2021 ZHANG Wei. All rights reserved.
 *
 * Author(s): ZHANG Wei <zhangwei@sustech.edu.cn>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include "constants.h"
#include "par_t.h"
// blk_t.h contains most other headers
#include "blk_t.h"

#include "media_discrete_model.h"
#include "drv_rk_curv_col.h"
#include "cuda_common.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

//-------------------------------------------------------------------------------
// initial gpu device before start MPI
//-------------------------------------------------------------------------------
  setDeviceBeforeInit();
//-------------------------------------------------------------------------------
// get commond-line argument
//-------------------------------------------------------------------------------

  // argc checking
  if (argc < 2) {
    fprintf(stdout,"usage: main_curv_col_el_2d <par_file> <opt: verbose>\n");
    exit(1);
  }

  //strncpy(par_fname, argv[1], sizeof(argv[1]));
  par_fname = argv[1];

  if (argc >= 3) {
    verbose = atoi(argv[2]); // verbose number
    fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
  }

  fprintf(stdout,"par file =  %s\n", par_fname); fflush(stdout);

  // read par

  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_read_from_file(par_fname, par, verbose);

  if (verbose>0) par_print(par);

//-------------------------------------------------------------------------------
// init blk_t
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"create blk ...\n"); 

  // malloc blk
  blk_t *blk = (blk_t *) malloc(sizeof(blk_t));

  // malloc inner vars
  blk_init(blk, verbose);

  fd_t            *fd            = blk->fd    ;
  gd_t            *gdcurv        = blk->gd;
  gdcurv_metric_t *gdcurv_metric = blk->gdcurv_metric;
  md_t            *md            = blk->md;
  wav_t           *wav           = blk->wav;
  src_t           *src           = blk->src;
  bdry_t          *bdry          = blk->bdry;
  iorecv_t        *iorecv        = blk->iorecv;
  ioline_t        *ioline        = blk->ioline;
  iosnap_t        *iosnap        = blk->iosnap;

  // set up fd_t
  //    not support selection scheme by par file yet
  if (verbose>0) fprintf(stdout,"set scheme ...\n"); 
  fd_set_macdrp(fd);

  // set gdinfo
  gd_indx_set(gdcurv, 
              par->number_of_total_grid_points_x,
              par->number_of_total_grid_points_z,
              fd->fdx_nghosts,
              fd->fdz_nghosts,
              verbose);

  // set str in blk
  blk_set_output(blk, 
                 par->output_dir,
                 par->grid_export_dir,
                 par->media_export_dir,
                 verbose);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"allocate grid vars ...\n"); 

  // malloc var in gdcurv
  gd_curv_init(gdcurv);

  // malloc var in gdcurv_metric
  gd_curv_metric_init(gdcurv, gdcurv_metric);

  // generate grid coord
  switch (par->grid_generation_itype)
  {
    case PAR_GRID_CARTESIAN : {

      fprintf(stdout,"generate cartesian grid in code ...\n"); 

      float dx = par->cartesian_grid_stepsize[0];
      float dz = par->cartesian_grid_stepsize[1];

      float x0 = par->cartesian_grid_origin[0];
      float z0 = par->cartesian_grid_origin[1];

      gd_curv_gen_cart(gdcurv,dx,x0,dz,z0);

      break;
    }
    case PAR_GRID_IMPORT : {

      fprintf(stdout,"import grid vars ...\n"); 
      gd_curv_coord_import(gdcurv, par->grid_import_dir);

      break;
    }
  }

  // cal min/max of this thread
  gd_curv_set_minmax(gdcurv);

  // generate topo over all the domain
  //ierr = gd_curv_topoall_generate();

  // output
  if (par->is_export_grid==1)
  {
    fprintf(stdout,"export coord to file ...\n"); 
    gd_curv_coord_export(gdcurv,
                         blk->grid_export_dir);
  } else {
    fprintf(stdout,"do not export coord\n"); 
  }
  fprintf(stdout, " --> done\n"); fflush(stdout);

  // cal metrics and output for QC
  switch (par->metric_method_itype)
  {
    case PAR_METRIC_CALCULATE : {

      if (verbose>0) fprintf(stdout,"calculate metrics ...\n"); 
      gd_curv_metric_cal(gdcurv,
                         gdcurv_metric,
                         fd->fdc_len,
                         fd->fdc_indx,
                         fd->fdc_coef);

      break;
    }
    case PAR_METRIC_IMPORT : {

      if (verbose>0) fprintf(stdout,"import metric file ...\n"); 
      gd_curv_metric_import(gdcurv_metric, par->grid_import_dir);

      break;
    }
  }
  if (verbose>0) { fprintf(stdout, " --> done\n"); fflush(stdout); }

  // export metric
  if (par->is_export_metric==1)
  {
    if (verbose>0) fprintf(stdout,"export metric to file ...\n"); 
    gd_curv_metric_export(gdcurv,gdcurv_metric,
                          blk->grid_export_dir);
  } else {
    if (verbose>0) fprintf(stdout,"do not export metric\n"); 
  }
  if (verbose>0) { fprintf(stdout, " --> done\n"); fflush(stdout); }

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (verbose>0) {fprintf(stdout,"allocate media vars ...\n"); fflush(stdout);}
  md_init(gdcurv, md, par->media_itype,  par->visco_itype);

  // read or discrete velocity model
  switch (par->media_input_itype)
  {
    case PAR_MEDIA_CODE : {

      if (verbose>0) fprintf(stdout,"generate simple medium in code ...\n"); 

      if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
        md_gen_test_ac_iso(md);
      }

      if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
        md_gen_test_el_iso(md);
      }

      if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI) {
        md_gen_test_el_vti(md);
      }

      if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
        md_gen_test_el_aniso(md);
      }

      if (md->visco_type == CONST_VISCO_GRAVES_QS) {
        md_gen_test_Qs(md, par->visco_Qs_freq);
      }

      break;
    }
    case PAR_MEDIA_IMPORT : {
      if (verbose>0) fprintf(stdout,"import discrete medium file ...\n"); 
      md_import(md, par->grid_import_dir);

      break;
    }

    case PAR_MEDIA_2LAY : {
      if (verbose>0) fprintf(stdout,"read and discretize layer medium file ...\n"); 

      if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
        media_layer2model_el_iso(md->lambda, md->mu, md->rho,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI) {
        media_layer2model_el_vti(md->rho,
                md->c11, md->c33, md->c55, md->c13,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
        media_layer2model_el_aniso(md->rho,
                md->c11, md->c13, md->c15, md->c33, md->c35, md->c55,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
        media_layer2model_ac_iso(md->rho, md->kappa,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      }

      break;
    }

    case PAR_MEDIA_2GRD : {
      if (verbose>0) fprintf(stdout,"read and descretize grid medium file ...\n"); 

      if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
        media_grid2model_el_iso(md->rho, md->lambda, md->mu,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                gdcurv->xmin, gdcurv->xmax, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI) {
        media_grid2model_el_vti(md->rho,
                md->c11, md->c33, md->c55, md->c13,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                gdcurv->xmin, gdcurv->xmax, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
        media_grid2model_el_aniso(md->rho,
                md->c11, md->c13, md->c15, md->c33, md->c35, md->c55,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                gdcurv->xmin, gdcurv->xmax, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      } else if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
        media_grid2model_ac_iso(md->rho, md->kappa,
                gdcurv->x2d, gdcurv->z2d,
                gdcurv->nx, gdcurv->nz, 
                gdcurv->xmin, gdcurv->xmax, 
                MEDIA_USE_CURV,
                par->media_input_file,
                par->equivalent_medium_method);
      }

      break;
    }
  }

  // export grid media
  if (par->is_export_media==1)
  {
    if (verbose>0) fprintf(stdout,"export discrete medium to file ...\n"); 

    md_export(gdcurv, md, blk->media_export_dir);
  } else {
    if (verbose>0) fprintf(stdout,"do not export medium\n"); 
  }

//-------------------------------------------------------------------------------
//-- estimate/check/set time step
//-------------------------------------------------------------------------------

  float   t0 = par->time_start;
  float   dt = par->size_of_time_step;
  int     nt_total = par->number_of_time_steps+1;

  if (par->time_check_stability==1)
  {
    float dtmax, dtmaxVp, dtmaxL;
    int   dtmaxi, dtmaxk;

    //-- estimate time step
    fprintf(stdout,"   estimate time step ...\n"); 
    blk_dt_esti_curv(gdcurv,md,fd->CFL,
            &dtmax, &dtmaxVp, &dtmaxL, &dtmaxi, &dtmaxk);
    
    //-- print for QC
    fprintf(stdout, "-> dtmax=%f, Vp=%f, L=%f, i=%d, k=%d\n",
            dtmax, dtmaxVp, dtmaxL, dtmaxi, dtmaxk);
    
    // check valid
    if (dtmax <= 0.0) {
       fprintf(stderr,"ERROR: maximum dt <= 0, stop running\n");
       exit(1);
    }

    //-- auto set stept
    if (dt < 0.0) {
       dt       = blk_keep_three_digi(dtmax);
       nt_total = (int) (par->time_window_length / dt + 0.5);

       fprintf(stdout, "-> Set dt       = %f according to maximum allowed value\n", dt);
       fprintf(stdout, "-> Set nt_total = %d\n", nt_total);
    }

    //-- if input dt, check value
    if (dtmax < dt) {
       fprintf(stdout, "Serious Error: dt=%f > dtmax=%f, stop!\n", dt, dtmax);
       exit(1);
    }
  }

//-------------------------------------------------------------------------------
//-- source import or locate on fly
//-------------------------------------------------------------------------------

  src_read_locate_file(gdcurv, src,
                       par->source_input_file,
                       t0, dt,
                       fd->num_rk_stages, fd->rk_rhs_time,
                       fd->fdx_max_half_len,
                       verbose);

  /*
  if (par->is_export_source==1)
  {
      ierr = src_export();
  }
  */

//-------------------------------------------------------------------------------
//-- allocate main var
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"allocate solver vars ...\n"); 
  if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO)
  {
    wav_ac_init(gdcurv, wav, fd->num_rk_stages);
  } else
  {
    wav_init(gdcurv, wav, fd->num_rk_stages);
  }

//-------------------------------------------------------------------------------
//-- setup output, may require coord info
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"setup output info ...\n"); 

  // receiver: need to do
  io_recv_read_locate(gdcurv, iorecv,
                      nt_total, wav->ncmp, par->in_station_file);

  // line
  io_line_locate(gdcurv, ioline,
                 wav->ncmp,
                 nt_total,
                 par->number_of_receiver_line,
                 par->receiver_line_index_start,
                 par->receiver_line_index_incre,
                 par->receiver_line_count,
                 par->receiver_line_name);
  
  // snapshot
  io_snapshot_locate(gdcurv, iosnap,
                     par->number_of_snapshot,
                     par->snapshot_name,
                     par->snapshot_index_start,
                     par->snapshot_index_count,
                     par->snapshot_index_incre,
                     par->snapshot_time_start,
                     par->snapshot_time_incre,
                     par->snapshot_save_velocity,
                     par->snapshot_save_stress,
                     par->snapshot_save_strain,
                     blk->output_dir);

//-------------------------------------------------------------------------------
//-- setup boundary
//-------------------------------------------------------------------------------

  if (verbose>0) fprintf(stdout,"setup boundary ...\n"); 

  bdry_init(bdry, gdcurv->nx, gdcurv->nz);

  //-- ade cfs-pml
  
  if (par->bdry_has_cfspml == 1)
  {
    if (verbose>0) fprintf(stdout,"setup ade cfs-pml ...\n"); 

    bdry_pml_set(gdcurv, wav, bdry,
                 par->cfspml_is_sides,
                 par->abs_num_of_layers,
                 par->cfspml_alpha_max,
                 par->cfspml_beta_max,
                 par->cfspml_velocity,
                 verbose);
  }

  //-- ablexp
  
  if (par->bdry_has_ablexp == 1)
  {
    if (verbose>0) fprintf(stdout,"setup sponge layer ...\n"); 

    bdry_ablexp_set(gdcurv, wav, bdry,
                    par->ablexp_is_sides,
                    par->abs_num_of_layers,
                    par->ablexp_velocity,
                    dt,
                    verbose);
  }

  //-- free surface preproc

  if (par->bdry_has_free == 1)
  {
    if (verbose>0) fprintf(stdout,"cal free surface matrix ...\n"); 

    bdry_free_set(gdcurv, bdry, par->free_is_sides, verbose);
  }

//-------------------------------------------------------------------------------
//-- qc
//-------------------------------------------------------------------------------
  
  fd_print(fd);

  blk_print(blk);

  gd_print(gdcurv);

  iosnap_print(iosnap);

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------
  
  // convert rho to 1 / rho to reduce number of arithmetic cal
  md_rho_to_slow(md->rho, md->siz_icmp);

  if (verbose>0) fprintf(stdout,"start solver ...\n"); 
  
  time_t t_start = time(NULL);
  
  drv_rk_curv_col_allstep(fd,gdcurv,gdcurv_metric,md,
                          src,bdry,
                          wav, 
                          iorecv,ioline,iosnap,
                          dt,nt_total,t0,
                          blk->output_dir,
                          par->check_nan_every_nummber_of_steps,
                          par->output_all,
                          verbose);
  
  time_t t_end = time(NULL);
  
  if (verbose>0) {
    fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));
  }

//-------------------------------------------------------------------------------
//-- save station and line seismo to sac
//-------------------------------------------------------------------------------

  io_recv_output_sac(iorecv,dt,wav->ncmp,wav->cmp_name,
                      src->evtnm,blk->output_dir,err_message);

  if(md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
    io_recv_output_sac_el_iso_strain(iorecv,md->lambda,md->mu,dt,
                          src->evtnm,blk->output_dir,err_message);
  }

  io_line_output_sac(ioline,dt,wav->cmp_name,src->evtnm,blk->output_dir);

  return 0;
}
