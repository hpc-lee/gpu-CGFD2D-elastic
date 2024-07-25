/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"

//
// malloc inner vars
//

int
blk_init(blk_t *blk, const int verbose)
{
  int ierr = 0;

  // alloc struct vars
  blk->fd            = (fd_t *)malloc(sizeof(fd_t));
  blk->gd            = (gd_t        *)malloc(sizeof(gd_t     ));
  blk->gdcurv_metric = (gdcurv_metric_t *)malloc(sizeof(gdcurv_metric_t));
  blk->md            = (md_t      *)malloc(sizeof(md_t     ));
  blk->wav           = (wav_t      *)malloc(sizeof(wav_t     ));
  blk->src           = (src_t      *)malloc(sizeof(src_t     ));
  blk->bdry          = (bdry_t     *)malloc(sizeof(bdry_t ));
  blk->iorecv        = (iorecv_t   *)malloc(sizeof(iorecv_t ));
  blk->ioline        = (ioline_t   *)malloc(sizeof(ioline_t ));
  blk->iosnap        = (iosnap_t   *)malloc(sizeof(iosnap_t ));

  sprintf(blk->name, "%s", "single");

  return ierr;
}

// set str
int
blk_set_output(blk_t *blk,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
               const int verbose)
{
  // set name
  //sprintf(blk->name, "%s", name);

  // output
  sprintf(blk->output_dir, "%s", output_dir);
  sprintf(blk->grid_export_dir, "%s", grid_export_dir);
  sprintf(blk->media_export_dir, "%s", media_export_dir);

  return 0;
}

int
blk_print(blk_t *blk)
{    
  int ierr = 0;

  fprintf(stdout, "\n-------------------------------------------------------\n");
  fprintf(stdout, "print blk %s:\n", blk->name);
  fprintf(stdout, "-------------------------------------------------------\n\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> ESTIMATE MEMORY INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "total memory size Byte: %20.5f  B\n", PSV->total_memory_size_Byte);
  //fprintf(stdout, "total memory size KB  : %20.5f KB\n", PSV->total_memory_size_KB  );
  //fprintf(stdout, "total memory size MB  : %20.5f MB\n", PSV->total_memory_size_MB  );
  //fprintf(stdout, "total memory size GB  : %20.5f GB\n", PSV->total_memory_size_GB  );
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> FOLDER AND FILE INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "   OutFolderName: %s\n", OutFolderName);
  //fprintf(stdout, "       EventName: %s\n", OutPrefix);
  //fprintf(stdout, "     LogFilename: %s\n", LogFilename);
  //fprintf(stdout, " StationFilename: %s\n", StationFilename);
  //fprintf(stdout, "  SourceFilename: %s\n", SourceFilename);
  //fprintf(stdout, "   MediaFilename: %s\n", MediaFilename);
  //fprintf(stdout, "\n");

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> media info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //if (blk->media_type == MEDIA_TYPE_LAYER)
  //{
  //    strcpy(str, "layer");
  //}
  //else if (blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    strcpy(str, "grid");
  //}
  //fprintf(stdout, " media_type = %s\n", str);
  //if(blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    fprintf(stdout, "\n --> the media filename is:\n");
  //    fprintf(stdout, " velp_file  = %s\n", blk->fnm_velp);
  //    fprintf(stdout, " vels_file  = %s\n", blk->fnm_vels);
  //    fprintf(stdout, " rho_file   = %s\n", blk->fnm_rho);
  //}
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> source info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, " number_of_force  = %d\n", blk->number_of_force);
  //if(blk->number_of_force > 0)
  //{
  //    fprintf(stdout, " force_source           x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_force; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->force_coord[indx], blk->force_coord[indx+1],
  //                blk->force_shift[indx], blk->force_shift[indx+1],
  //                blk->force_indx [indx], blk->force_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "\n");
  //fprintf(stdout, " number_of_moment = %d\n", blk->number_of_moment);
  //if(blk->number_of_moment > 0)
  //{
  //    fprintf(stdout, " moment_source          x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_moment; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->moment_coord[indx], blk->moment_coord[indx+1],
  //                blk->moment_shift[indx], blk->moment_shift[indx+1],
  //                blk->moment_indx [indx], blk->moment_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> boundary layer information:\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //ierr = boundary_id2type(type1, blk->boundary_type[0], errorMsg);
  //ierr = boundary_id2type(type2, blk->boundary_type[1], errorMsg);
  //ierr = boundary_id2type(type3, blk->boundary_type[2], errorMsg);
  //ierr = boundary_id2type(type4, blk->boundary_type[3], errorMsg);
  //fprintf(stdout, " boundary_type         = %10s%10s%10s%10s\n", 
  //        type1, type2, type3, type4);
  //fprintf(stdout, " boundary_layer_number = %10d%10d%10d%10d\n", 
  //        blk->boundary_layer_number[0], blk->boundary_layer_number[1], 
  //        blk->boundary_layer_number[2], blk->boundary_layer_number[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " absorb_velocity       = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->absorb_velocity[0], blk->absorb_velocity[1], blk->absorb_velocity[2], 
  //        blk->absorb_velocity[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " CFS_alpha_max         = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_alpha_max[0], blk->CFS_alpha_max[1], blk->CFS_alpha_max[2], 
  //        blk->CFS_alpha_max[3]);
  //fprintf(stdout, " CFS_beta_max          = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_beta_max[0], blk->CFS_beta_max[1], blk->CFS_beta_max[2], 
  //        blk->CFS_beta_max[3]);
  
  return ierr;
}

/*********************************************************************
 * estimate dt
 *********************************************************************/

int
blk_dt_esti_curv(gd_t *gdcurv, md_t *md,
                 float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
                 int *dtmaxi, int *dtmaxk)
{
  int ierr = 0;

  float dtmax_local = 1.0e10;
  float Vp;

  float * x2d = gdcurv->x2d;
  float * z2d = gdcurv->z2d;

  for (int k = gdcurv->nk1; k <= gdcurv->nk2; k++)
  {
      for (int i = gdcurv->ni1; i <= gdcurv->ni2; i++)
      {
        size_t iptr = i + k * gdcurv->siz_iz;

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( (md->lambda[iptr] + 2.0 * md->mu[iptr]) / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI) {
          float Vpv = sqrt( md->c33[iptr] / md->rho[iptr] );
          float Vph = sqrt( md->c11[iptr] / md->rho[iptr] );
          Vp = Vph > Vpv ? Vph : Vpv;
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
          // need to implement accurate solution
          Vp = sqrt( md->c11[iptr] / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
          Vp = sqrt( md->kappa[iptr] / md->rho[iptr] );
        }

        float dtLe = 1.0e20;
        float x0 = x2d[iptr];
        float z0 = z2d[iptr];

        // min L to 8 adjacent planes
        for (int kk = -1; kk <=1; kk++) {
            for (int ii = -1; ii <= 1; ii++) {
              if (ii != 0 && kk != 0)
              {
                float p1[] = { x2d[iptr-ii], z2d[iptr-ii] };
                float p2[] = { x2d[iptr-kk*gdcurv->siz_iz],
                               z2d[iptr-kk*gdcurv->siz_iz] };

                float L = fdlib_math_dist_point2line(x0, z0, p1, p2);

                if (dtLe > L) dtLe = L;
              }
          }
        }

        // convert to dt
        float dt_point = CFL / Vp * dtLe;

        // if smaller
        if (dt_point < dtmax_local) {
          dtmax_local = dt_point;
          *dtmaxi = i;
          *dtmaxk = k;
          *dtmaxVp = Vp;
          *dtmaxL  = dtLe;
        }
      } // i
  } //k

  *dtmax = dtmax_local;

  return ierr;
}

int
blk_dt_esti_cart(gd_t *gdcart, md_t *md,
                 float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
                 int *dtmaxi, int *dtmaxk)
{
  int ierr = 0;

  float dtmax_local = 1.0e10;
  float Vp;

  float dx = gdcart->dx;
  float dz = gdcart->dz;

  // length to plane
  float p1[] = {  dx, 0.0 };
  float p2[] = { 0.0,  dz };

  float dtLe = fdlib_math_dist_point2line(0.0,0.0, p1, p2);

  for (int k = gdcart->nk1; k <= gdcart->nk2; k++)
  {
      for (int i = gdcart->ni1; i <= gdcart->ni2; i++)
      {
        size_t iptr = i + k * gdcart->siz_iz;

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( (md->lambda[iptr] + 2.0 * md->mu[iptr]) / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_VTI) {
          float Vpv = sqrt( md->c33[iptr] / md->rho[iptr] );
          float Vph = sqrt( md->c11[iptr] / md->rho[iptr] );
          Vp = Vph > Vpv ? Vph : Vpv;
        } else if (md->medium_type == CONST_MEDIUM_ACOUSTIC_ISO) {
          Vp = sqrt( md->kappa[iptr] / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
          Vp = sqrt( (md->c11[iptr]) / md->rho[iptr] );
        } else {
          fprintf(stderr,"ERROR: medium type is not implemented\n");
          exit(1);
        }

        // convert to dt
        float dt_point = CFL / Vp * dtLe;

        // if smaller
        if (dt_point < dtmax_local) {
          dtmax_local = dt_point;
          *dtmaxi = i;
          *dtmaxk = k;
          *dtmaxVp = Vp;
        }

      } // i
  } //k

  *dtmax  = dtmax_local;
  *dtmaxL = dtLe;

  return ierr;
}

float
blk_keep_three_digi(float dt)
{
  char str[40];
  float dt_2;

  sprintf(str, "%9.7e", dt);

  for (int i = 3; i < 9; i++)
  {
    str[i] = '0';
  }

  sscanf(str, "%f", &dt_2);
  
  return dt_2;
}
