#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"
#include "sacLib.h"

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "constants.h"
#include "fd_t.h"
#include "io_funcs.h"
#include "cuda_common.h"

//#define M_NCERR(ierr) {fprintf(stderr,"sv_ nc error: %s\n", nc_strerror(ierr)); exit(1);}
#ifndef M_NCERR
#define M_NCERR {fprintf(stderr,"io nc error\n"); exit(1);}
#endif

/*
 * read in station list file and locate station
 */
int
io_recv_read_locate(gd_t *gd,
                    iorecv_t  *iorecv,
                    int       nt_total,
                    int       num_of_vars,
                    char *in_filenm)
{
  FILE *fp;
  char line[500];

  if (!(fp = fopen (in_filenm, "rt")))
	{
	    fprintf (stdout, "Cannot open input file %s\n", in_filenm);
	    fflush (stdout);
	    return 1;
	}

  // number of station
  int nr;

  io_get_nextline(fp, line, 500);
  sscanf(line, "%d", &nr);

  // check fail in the future
  iorecv_one_t *recvone = (iorecv_one_t *)malloc(nr * sizeof(iorecv_one_t));

  // read coord and locate

  int ir=0;
  int nr_this = 0; // in this thread
  int is_coord;
  int is_depth;

  for (ir=0; ir<nr; ir++)
  {
    float rx, rz;
    float rx_inc, rz_inc; // shift in computational space grid
    int ix, iz; // local index

    // read one line
    io_get_nextline(fp, line, 500);

    // get values
    sscanf(line, "%s %d %d %g %g", 
              recvone[ir].name, &is_coord, &is_depth, &rx, &rz);
    
    // by grid index
    if (is_coord == 0)
    {
      // add ghosts, to local point
      rx = rx + gd->fdx_nghosts;
      rz = rz + gd->fdz_nghosts;
      // if sz is relative to surface, convert to normal index
      if (is_depth == 1) {
        rz = gd->nk2 - rz + gd->fdz_nghosts;
      }

      // do not take nearest value, but use smaller value
      ix = (int) (rx + 0.0);
      iz = (int) (rz + 0.0);
      rx_inc = rx - ix;
      rz_inc = rz - iz;
    }
    // by axis
    else
    {
      // convert coord to global index
      if (gd->type == GD_TYPE_CURV)
      {
        // if rz is depth, convert to axis when it is in this thread
        if (is_depth == 1) {
          gd_curv_depth_to_axis(gd,rx,&rz);
        }
        gd_curv_coord_to_local_indx(gd,rx,rz,
                               &ix,&iz,&rx_inc,&rz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        // if sz is depth, convert to axis
        if (is_depth == 1) {
          rz = gd->z1d[gd->nk2] - rz;
        }
        gd_cart_coord_to_local_indx(gd,rx,rz,
                               &ix,&iz,&rx_inc,&rz_inc);
      }

      if (rx_inc < 0.0) {
        rx_inc = 1.0 + rx_inc;
        ix -= 1;
      }
      if (rz_inc < 0.0) {
        rz_inc = 1.0 + rz_inc;
        iz -= 1;
      }
    }

    if (gd_lindx_is_inner(ix,iz,gd) == 1)
    {
      // index, to get coord
      if (is_coord == 0)
      {
        rx = gd_coord_get_x(gd,ix,iz);
        rz = gd_coord_get_z(gd,ix,iz);
      }

      iorecv_one_t *this_recv = recvone + nr_this;

      sprintf(this_recv->name, "%s", recvone[ir].name);

      // get coord
      this_recv->x = rx;
      this_recv->z = rz;
      // set point and shift
      this_recv->i=ix;
      this_recv->k=iz;
      this_recv->di = rx_inc;
      this_recv->dk = rz_inc;

      this_recv->indx1d[0] = ix   + iz     * gd->siz_iz;
      this_recv->indx1d[1] = ix+1 + iz     * gd->siz_iz;
      this_recv->indx1d[2] = ix   + (iz+1) * gd->siz_iz;
      this_recv->indx1d[3] = ix+1 + (iz+1) * gd->siz_iz;

      nr_this += 1;
    }
  }

  fclose(fp);

  iorecv->total_number = nr_this;
  iorecv->recvone      = recvone;
  iorecv->max_nt       = nt_total;
  iorecv->ncmp         = num_of_vars;

  // malloc seismo
  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    recvone = iorecv->recvone + ir;
    recvone->seismo = (float *) malloc(num_of_vars * nt_total * sizeof(float));
  }

  return 0;
}

int io_line_locate(gd_t *gd,
                   ioline_t *ioline,
                   int    num_of_vars,
                   int    nt_total,
                   int    number_of_receiver_line,
                   int   *receiver_line_index_start,
                   int   *receiver_line_index_incre,
                   int   *receiver_line_count,
                   char **receiver_line_name)
{
  int ierr = 0;

  // init
  ioline->num_of_lines  = 0;
  ioline->max_nt        = nt_total;
  ioline->ncmp          = num_of_vars;

  // alloc as max num to keep nr and seq values, easy for second round
  ioline->line_nr  = (int *) malloc(number_of_receiver_line * sizeof(int));
  ioline->line_seq = (int *) malloc(number_of_receiver_line * sizeof(int));

  // first run to count line and nr
  for (int n=0; n < number_of_receiver_line; n++)
  {
    int nr = 0;
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*CONST_NDIM+0] 
                 + ipt * receiver_line_index_incre[n*CONST_NDIM];

      int gk = receiver_line_index_start[n*CONST_NDIM+1] 
                 + ipt * receiver_line_index_incre[n*CONST_NDIM+1];

      if (gd_pindx_is_inner(gi,gk,gd) == 1)
      {
        nr += 1;
      }
    }

    // if any receiver of this line in this thread
    if (nr>0)
    {
      ioline->line_nr [ ioline->num_of_lines ] = nr;
      ioline->line_seq[ ioline->num_of_lines ] = n;
      ioline->num_of_lines += 1;
    }
  }

  // alloc
  if (ioline->num_of_lines>0)
  {
    ioline->line_name   = (char **)fdlib_mem_malloc_2l_char(ioline->num_of_lines,
                                    CONST_MAX_STRLEN, "io_line_locate");

    ioline->recv_seq    = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    //ioline->recv_indx   = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    ioline->recv_iptr   = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    ioline->recv_x  = (float **) malloc(ioline->num_of_lines * sizeof(float*));
    ioline->recv_z  = (float **) malloc(ioline->num_of_lines * sizeof(float*));
    ioline->recv_seismo = (float **) malloc(ioline->num_of_lines * sizeof(float*));

    for (int n=0; n < ioline->num_of_lines; n++)
    {
      int nr = ioline->line_nr[n];
      //ioline->recv_indx[n] = (int *)malloc(nr * CONST_NDIM * sizeof(int)); 
      ioline->recv_seq [n]  = (int *)malloc( nr * sizeof(int) ); 
      ioline->recv_iptr[n]  = (int *)malloc( nr * sizeof(int) ); 
      ioline->recv_x[n] = (float *)malloc( nr * sizeof(float) );
      ioline->recv_z[n] = (float *)malloc( nr * sizeof(float) );
      ioline->recv_seismo[n] = (float *)malloc(
                                nr * num_of_vars * nt_total * sizeof(float) );
    }
  }

  // second run for value
  //  only loop lines in this thread
  for (int m=0; m < ioline->num_of_lines; m++)
  {
    int n = ioline->line_seq[m];

    sprintf(ioline->line_name[m], "%s", receiver_line_name[n]);

    int ir = 0;
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*CONST_NDIM+0] + ipt * receiver_line_index_incre[n*CONST_NDIM  ];
      int gk = receiver_line_index_start[n*CONST_NDIM+1] + ipt * receiver_line_index_incre[n*CONST_NDIM+1];

      if (gd_pindx_is_inner(gi,gk,gd) == 1)
      {
        int i = gi + gd->fdx_nghosts;
        int k = gk + gd->fdz_nghosts;

        int iptr = i + k * gd->siz_iz;

        ioline->recv_seq [m][ir] = ipt;
        ioline->recv_iptr[m][ir] = iptr;

        ioline->recv_x[m][ir] = gd_coord_get_x(gd,i,k);
        ioline->recv_z[m][ir] = gd_coord_get_z(gd,i,k);

        ir += 1;
      }
    }
  }

  return ierr;
}

int
io_snapshot_locate(gd_t *gd,
                   iosnap_t *iosnap,
                    int  number_of_snapshot,
                    char **snapshot_name,
                    int *snapshot_index_start,
                    int *snapshot_index_count,
                    int *snapshot_index_incre,
                    int *snapshot_time_start,
                    int *snapshot_time_incre,
                    int *snapshot_save_velocity,
                    int *snapshot_save_stress,
                    int *snapshot_save_strain,
                    char *output_dir)
{
  // malloc to max, num of snap will not be large
  if (number_of_snapshot > 0)
  {
    iosnap->fname = (char **) fdlib_mem_malloc_2l_char(number_of_snapshot,
                                    CONST_MAX_STRLEN,"snap_fname");
    iosnap->i1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->k1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->ni = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->nk = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->di = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->dk = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->it1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->dit = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_vel    = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_stress = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_strain = (int *) malloc(number_of_snapshot * sizeof(int));

    iosnap->i1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->k1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
  }

  // init

  iosnap->siz_max_wrk = 0;

  int isnap = 0;

  for (int n=0; n < number_of_snapshot; n++)
  {
    int iptr0 = n*CONST_NDIM;

    // scan output k-index in this proc
    int gk1 = -1; int ngk =  0; int k_in_nc = 0;
    for (int n3=0; n3<snapshot_index_count[iptr0+1]; n3++)
    {
      int gk = snapshot_index_start[iptr0+1] + n3 * snapshot_index_incre[iptr0+1];
      if (gd_pindx_is_inner_k(gk,gd) == 1)
      {
        // first valid k
        if (gk1 == -1) {
          gk1 = gk;
          k_in_nc = n3;
        }
        ngk++;
      }
      if (gk > gd->nk-1) break; // no need to larger k
    }

    // scan output i-index in this proc
    int gi1 = -1; int ngi =  0; int i_in_nc = 0;
    for (int n1=0; n1<snapshot_index_count[iptr0+0]; n1++)
    {
      int gi = snapshot_index_start[iptr0+0] + n1 * snapshot_index_incre[iptr0+0];
      if (gd_pindx_is_inner_i(gi,gd) == 1)
      {
        if (gi1 == -1) {
          gi1 = gi;
          i_in_nc = n1;
        }
        ngi++;
      }
      if (gi > gd->ni-1) break;
    }

    // if in this proc
    if (ngi>0 && ngk>0)
    {
      iosnap->i1[isnap]  = gi1 + gd->fdx_nghosts;
      iosnap->k1[isnap]  = gk1 + gd->fdz_nghosts;
      iosnap->ni[isnap]  = ngi;
      iosnap->nk[isnap]  = ngk;
      iosnap->di[isnap]  = snapshot_index_incre[iptr0+0];
      iosnap->dk[isnap]  = snapshot_index_incre[iptr0+1];

      iosnap->it1[isnap]  = snapshot_time_start[n];
      iosnap->dit[isnap]  = snapshot_time_incre[n];

      iosnap->out_vel   [isnap] = snapshot_save_velocity[n];
      iosnap->out_stress[isnap] = snapshot_save_stress[n];
      iosnap->out_strain[isnap] = snapshot_save_strain[n];

      iosnap->i1_to_glob[isnap] = i_in_nc;
      iosnap->k1_to_glob[isnap] = k_in_nc;

      sprintf(iosnap->fname[isnap],"%s/%s.nc",output_dir,
                                                 snapshot_name[n]);

      // for max wrk
      size_t snap_siz =  ngi * ngk;
      iosnap->siz_max_wrk = snap_siz > iosnap->siz_max_wrk ? 
                            snap_siz : iosnap->siz_max_wrk;

      isnap += 1;
    } // if in this
  } // loop all snap

  iosnap->num_of_snap = isnap;

  return 0;
}

int
io_snap_nc_create(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc)
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  char **snap_fname = iosnap->fname;

  iosnap_nc->num_of_snap = num_of_snap;
  iosnap_nc->ncid = (int *)malloc(num_of_snap*sizeof(int));
  iosnap_nc->timeid = (int *)malloc(num_of_snap*sizeof(int));

  iosnap_nc->varid_V = (int *)malloc(num_of_snap*CONST_NDIM*sizeof(int));
  iosnap_nc->varid_T = (int *)malloc(num_of_snap*3 *sizeof(int));
  iosnap_nc->varid_E = (int *)malloc(num_of_snap*3 *sizeof(int));

  // will be used in put step
  iosnap_nc->cur_it = (int *)malloc(num_of_snap*sizeof(int));
  for (int n=0; n<num_of_snap; n++) {
    iosnap_nc->cur_it[n] = 0;
  }

  int *ncid   = iosnap_nc->ncid;
  int *timeid = iosnap_nc->timeid;
  int *varid_V = iosnap_nc->varid_V;
  int *varid_T = iosnap_nc->varid_T;
  int *varid_E = iosnap_nc->varid_E;

  for (int n=0; n<num_of_snap; n++)
  {
    int dimid[CONST_NDIM+1];
    int snap_i1  = iosnap->i1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dk  = iosnap->dk[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];

    ierr = nc_create(snap_fname[n], NC_CLOBBER, &ncid[n]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "time", NC_UNLIMITED, &dimid[0]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "k", snap_nk     , &dimid[1]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "i", snap_ni     , &dimid[2]);
    handle_nc_err(ierr);
    // time var
    ierr = nc_def_var(ncid[n], "time", NC_FLOAT, 1, dimid+0, &timeid[n]);
    handle_nc_err(ierr);
    // other vars
    if (snap_out_V==1) {
       ierr = nc_def_var(ncid[n],"Vx",NC_FLOAT,CONST_NDIM+1,dimid,&varid_V[n*CONST_NDIM+0]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Vz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_V[n*CONST_NDIM+1]);
       handle_nc_err(ierr);
    }
    if (snap_out_T==1) {
       ierr = nc_def_var(ncid[n],"Txx",NC_FLOAT,CONST_NDIM+1,dimid,&varid_T[n*3+0]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Tzz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_T[n*3+1]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Txz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_T[n*3+2]);
       handle_nc_err(ierr);
    }
    if (snap_out_E==1) {
       ierr = nc_def_var(ncid[n],"Exx",NC_FLOAT,CONST_NDIM+1,dimid,&varid_E[n*3 +0]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Ezz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_E[n*3 +1]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Exz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_E[n*3 +2]);
       handle_nc_err(ierr);
    }
    // attribute: index in output snapshot, index w ghost in thread
    int g_start[] = { iosnap->i1_to_glob[n],
                      iosnap->k1_to_glob[n] };
    nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_to_snapshot_output",
                   NC_INT,CONST_NDIM,g_start);

    int l_start[] = { snap_i1, snap_k1 };
    nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_in_this_thread_with_ghosts",
                   NC_INT,CONST_NDIM,l_start);

    int l_count[] = { snap_di, snap_dk };
    nc_put_att_int(ncid[n],NC_GLOBAL,"index_stride_in_this_thread",
                   NC_INT,CONST_NDIM,l_count);

    ierr = nc_enddef(ncid[n]);
    handle_nc_err(ierr);
  } // loop snap

  return ierr;
}

/*
 * 
 */

int
io_snap_nc_put(iosnap_t *iosnap,
               iosnap_nc_t *iosnap_nc,
               gd_t    *gd,
               md_t    *md,
               wav_t   *wav,
               float * w_end_d,
               float * buff,
               int   nt_total,
               int   it,
               float time,
               int is_run_out_vel,     // for stg, out vel and stress at sep call
               int is_run_out_stress,  // 
               int is_incr_cur_it)     // for stg, should output cur_it once
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  size_t siz_iz = gd->siz_iz;
  size_t siz_icmp = gd->siz_icmp;

  for (int n=0; n<num_of_snap; n++)
  {
    int snap_i1  = iosnap->i1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dk  = iosnap->dk[n];

    int snap_it1 = iosnap->it1[n];
    int snap_dit = iosnap->dit[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];

    int snap_it_mod = (it - snap_it1) % snap_dit;
    int snap_it_num = (it - snap_it1) / snap_dit;
    int snap_nt_total = (nt_total - snap_it1) / snap_dit;

    int snap_max_num = snap_ni * snap_nk;
    float *buff_d = NULL;

    if (it>=snap_it1 && snap_it_num<=snap_nt_total && snap_it_mod==0)
    {
      size_t startp[] = { iosnap_nc->cur_it[n], 0, 0 };
      size_t countp[] = { 1, snap_nk,  snap_ni };
      size_t start_tdim = iosnap_nc->cur_it[n];

      // put time var
      nc_put_var1_float(iosnap_nc->ncid[n],iosnap_nc->timeid[n],&start_tdim,&time);
      int size = sizeof(float)*snap_max_num;
      buff_d = (float *) cuda_malloc(size);
      dim3 block(8,8);
      dim3 grid;
      grid.x = (snap_ni+block.x-1)/block.x;
      grid.y = (snap_nk+block.y-1)/block.y;

      // vel
      if (is_run_out_vel == 1 && snap_out_V==1)
      {
        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Vx_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+0*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+0],
              startp,countp,buff+0*siz_icmp);

        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Vz_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+1*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+1],
              startp,countp,buff+1*siz_icmp);
      }

      if (is_run_out_stress==1 && snap_out_T==1)
      {
        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Txx_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+2*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*3+0],
              startp,countp,buff+2*siz_icmp);

        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Tzz_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+3*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*3+1],
              startp,countp,buff+3*siz_icmp);

        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Txz_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+4*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*3+2],
              startp,countp,buff+4*siz_icmp);
      }

      if (is_run_out_stress==1 && snap_out_E==1 &&
          md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
      {
        // if snap_out_T==0, output T to calculate E
        if (is_run_out_stress==1 && snap_out_T==0)
        {
          io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Txx_pos,
                   siz_iz,snap_i1,snap_ni,snap_di,
                   snap_k1,snap_nk,snap_dk,buff_d);
          CUDACHECK(cudaMemcpy(buff+2*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));

          io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Tzz_pos,
                   siz_iz,snap_i1,snap_ni,snap_di,
                   snap_k1,snap_nk,snap_dk,buff_d);
          CUDACHECK(cudaMemcpy(buff+3*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));

          io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Txz_pos,
                   siz_iz,snap_i1,snap_ni,snap_di,
                   snap_k1,snap_nk,snap_dk,buff_d);
          CUDACHECK(cudaMemcpy(buff+4*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        }
        // convert to strain
        io_snap_stress_to_strain_eliso(md->lambda,md->mu,
                                       buff + 2*siz_icmp,   //Txx
                                       buff + 3*siz_icmp,   //Tzz
                                       buff + 4*siz_icmp,   //Txz
                                       buff + 5*siz_icmp,   //Exx
                                       buff + 6*siz_icmp,   //Ezz
                                       buff + 7*siz_icmp,   //Exz
                                       siz_iz,snap_i1,snap_ni,snap_di,
                                       snap_k1,snap_nk,snap_dk);
        // export
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*3+0],
              startp,countp,buff + 5*siz_icmp);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*3+1],
              startp,countp,buff + 6*siz_icmp);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*3+2],
              startp,countp,buff + 7*siz_icmp);

      }

      if (is_incr_cur_it == 1) {
        iosnap_nc->cur_it[n] += 1;
      }

      CUDACHECK(cudaFree(buff_d));

    } // if it
  } // loop snap

  return ierr;
}

int
io_snap_nc_create_ac(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc)
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  char **snap_fname = iosnap->fname;

  iosnap_nc->num_of_snap = num_of_snap;
  iosnap_nc->ncid = (int *)malloc(num_of_snap*sizeof(int));
  iosnap_nc->timeid = (int *)malloc(num_of_snap*sizeof(int));

  iosnap_nc->varid_V = (int *)malloc(num_of_snap*CONST_NDIM*sizeof(int));
  iosnap_nc->varid_T = (int *)malloc(num_of_snap*1*sizeof(int));
  iosnap_nc->varid_E = (int *)malloc(num_of_snap*1*sizeof(int));

  // will be used in put step
  iosnap_nc->cur_it = (int *)malloc(num_of_snap*sizeof(int));
  for (int n=0; n<num_of_snap; n++) {
    iosnap_nc->cur_it[n] = 0;
  }

  int *ncid   = iosnap_nc->ncid;
  int *timeid = iosnap_nc->timeid;
  int *varid_V = iosnap_nc->varid_V;
  int *varid_T = iosnap_nc->varid_T;
  int *varid_E = iosnap_nc->varid_E;

  for (int n=0; n<num_of_snap; n++)
  {
    int dimid[CONST_NDIM+1];
    int snap_i1  = iosnap->i1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dk  = iosnap->dk[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];

    ierr = nc_create(snap_fname[n], NC_CLOBBER, &ncid[n]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "time", NC_UNLIMITED, &dimid[0]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "k", snap_nk     , &dimid[1]);
    handle_nc_err(ierr);
    ierr = nc_def_dim(ncid[n], "i", snap_ni     , &dimid[2]);
    handle_nc_err(ierr);
    // time var
    ierr = nc_def_var(ncid[n], "time", NC_FLOAT, 1, dimid+0, &timeid[n]);
    handle_nc_err(ierr);
    // other vars
    if (snap_out_V==1) {
       ierr = nc_def_var(ncid[n],"Vx",NC_FLOAT,CONST_NDIM+1,dimid,&varid_V[n*CONST_NDIM+0]);
       handle_nc_err(ierr);
       ierr = nc_def_var(ncid[n],"Vz",NC_FLOAT,CONST_NDIM+1,dimid,&varid_V[n*CONST_NDIM+1]);
       handle_nc_err(ierr);
    }
    if (snap_out_T==1) {
       ierr = nc_def_var(ncid[n],"Tii",NC_FLOAT,CONST_NDIM+1,dimid,&varid_T[n]);
       handle_nc_err(ierr);
    }
    if (snap_out_E==1) {
       ierr = nc_def_var(ncid[n],"Eii",NC_FLOAT,CONST_NDIM+1,dimid,&varid_E[n]);
       handle_nc_err(ierr);
    }
    // attribute: index in output snapshot, index w ghost in thread
    int g_start[] = { iosnap->i1_to_glob[n],
                      iosnap->k1_to_glob[n] };
    nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_to_snapshot_output",
                   NC_INT,CONST_NDIM,g_start);

    int l_start[] = { snap_i1, snap_k1 };
    nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_in_this_thread_with_ghosts",
                   NC_INT,CONST_NDIM,l_start);

    int l_count[] = { snap_di, snap_dk };
    nc_put_att_int(ncid[n],NC_GLOBAL,"index_stride_in_this_thread",
                   NC_INT,CONST_NDIM,l_count);

    ierr = nc_enddef(ncid[n]);
    handle_nc_err(ierr);
  } // loop snap

  return ierr;
}

int
io_snap_nc_put_ac(iosnap_t *iosnap,
               iosnap_nc_t *iosnap_nc,
               gd_t    *gd,
               md_t    *md,
               wav_t   *wav,
               float * w_end_d,
               float * buff,
               int   nt_total,
               int   it,
               float time,
               int is_run_out_vel,     // for stg, out vel and stress at sep call
               int is_run_out_stress,  // 
               int is_incr_cur_it)     // for stg, should output cur_it once
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  size_t siz_iz   = gd->siz_iz;
  size_t siz_icmp = gd->siz_icmp;

  for (int n=0; n<num_of_snap; n++)
  {
    int snap_i1  = iosnap->i1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dk  = iosnap->dk[n];

    int snap_it1 = iosnap->it1[n];
    int snap_dit = iosnap->dit[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];

    int snap_it_mod = (it - snap_it1) % snap_dit;
    int snap_it_num = (it - snap_it1) / snap_dit;
    int snap_nt_total = (nt_total - snap_it1) / snap_dit;

    int snap_max_num = snap_ni * snap_nk;
    float *buff_d = NULL;

    if (it>=snap_it1 && snap_it_num<=snap_nt_total && snap_it_mod==0)
    {
      size_t startp[] = { iosnap_nc->cur_it[n], 0, 0 };
      size_t countp[] = { 1, snap_nk, snap_ni };
      size_t start_tdim = iosnap_nc->cur_it[n];

      // put time var
      nc_put_var1_float(iosnap_nc->ncid[n],iosnap_nc->timeid[n],&start_tdim,&time);
      int size = sizeof(float)*snap_max_num;
      buff_d = (float *) cuda_malloc(size);
      dim3 block(8,8);
      dim3 grid;
      grid.x = (snap_ni+block.x-1)/block.x;
      grid.y = (snap_nk+block.y-1)/block.y;

      // vel
      if (is_run_out_vel == 1 && snap_out_V==1)
      {
        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Vx_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+0*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+0],
              startp,countp,buff+0*siz_icmp);

        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Vz_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+1*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+1],
              startp,countp,buff+1*siz_icmp);
      }
      if (is_run_out_stress==1 && snap_out_T==1)
      {
        io_snap_pack_buff<<<grid, block>>> (w_end_d + wav->Txx_pos,
                 siz_iz,snap_i1,snap_ni,snap_di,
                 snap_k1,snap_nk,snap_dk,buff_d);
        CUDACHECK(cudaMemcpy(buff+2*siz_icmp,buff_d,size,cudaMemcpyDeviceToHost));
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n],
              startp,countp,buff+2*siz_icmp);
      }
      if (is_run_out_stress==1 && snap_out_E==1)
      {
        // need to implement
      }

      if (is_incr_cur_it == 1) {
        iosnap_nc->cur_it[n] += 1;
      }

    } // if it
  } // loop snap

  return ierr;
}

__global__ void
io_snap_pack_buff(float * var,
                  size_t siz_iz,
                  int starti,
                  int counti,
                  int increi,
                  int startk,
                  int countk,
                  int increk,
                  float * buff_d)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  size_t iz = blockIdx.y * blockDim.y + threadIdx.y;
  if(ix<counti && iz<countk)
  {
    size_t iptr_snap = ix + iz * counti;
    size_t i = starti + ix * increi;
    size_t k = startk + iz * increk;
    size_t iptr = i + k * siz_iz;
    buff_d[iptr_snap] =  var[iptr];
  }
  return;
}

int
io_snap_nc_close(iosnap_nc_t *iosnap_nc)
{
  for (int n=0; n < iosnap_nc->num_of_snap; n++)
  {
    nc_close(iosnap_nc->ncid[n]);
  }
  return 0;
}

int
io_snap_stress_to_strain_eliso(float *lam3d,
                               float *mu3d,
                               float *Txx,
                               float *Tzz,
                               float *Txz,
                               float *Exx,
                               float *Ezz,
                               float *Exz,
                               size_t siz_iz,
                               int starti,
                               int counti,
                               int increi,
                               int startk,
                               int countk,
                               int increk)
{
  size_t iptr_snap=0;
  size_t i,k,iptr,iptr_k;
  float lam,mu,E1,E2,E3,E0;

  for (int n3=0; n3<countk; n3++)
  {
    k = startk + n3 * increk;
    iptr_k = k * siz_iz;
    for (int n1=0; n1<counti; n1++)
    {
      i = starti + n1 * increi;
      iptr = i + iptr_k;
      iptr_snap = n1 + n3 * counti;

      lam = lam3d[iptr];
      mu  =  mu3d[iptr];
      
      E1 = (lam + mu) / (mu * ( 3.0 * lam + 2.0 * mu));
      E2 = - lam / ( 2.0 * mu * (3.0 * lam + 2.0 * mu));
      E3 = 1.0 / mu;

      E0 = E2 * (Txx[iptr_snap] + Tzz[iptr_snap]);

      Exx[iptr_snap] = E0 - (E2 - E1) * Txx[iptr_snap];
      Ezz[iptr_snap] = E0 - (E2 - E1) * Tzz[iptr_snap];
      Exz[iptr_snap] = 0.5 * E3 * Txz[iptr_snap];
    } //i
  } //k

  return 0;
}

int
io_recv_keep(iorecv_t *iorecv, float * w_end_d,
             float * buff, int it, int ncmp, int siz_icmp)
{
  float Lx1, Lx2, Lz1, Lz2;
  //CONST_NDIM_2 = 4, use 4 points interp
  int size = sizeof(float)*ncmp*CONST_NDIM_2;
  float *buff_d = (float *) cuda_malloc(size);
  size_t *indx1d_d = (size_t *) cuda_malloc(sizeof(size_t)*CONST_NDIM_2);
  dim3 block(32);
  dim3 grid;
  grid.x = (ncmp+block.x-1)/block.x;

  for (int n=0; n < iorecv->total_number; n++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + n;
    size_t *indx1d = this_recv->indx1d;
    CUDACHECK(cudaMemcpy(indx1d_d,indx1d,sizeof(size_t)*CONST_NDIM_2,cudaMemcpyHostToDevice));

    // get coef of linear interp
    Lx2 = this_recv->di; Lx1 = 1.0 - Lx2;
    Lz2 = this_recv->dk; Lz1 = 1.0 - Lz2;

    io_recv_line_interp_pack_buff<<<grid, block>>> (w_end_d, buff_d, ncmp, siz_icmp, indx1d_d);
    CUDACHECK(cudaMemcpy(buff,buff_d,size,cudaMemcpyDeviceToHost));

    for (int icmp=0; icmp < ncmp; icmp++)
    {
      int iptr_sta = icmp * iorecv->max_nt + it;
      this_recv->seismo[iptr_sta] =  buff[CONST_NDIM_2*icmp + 0] * Lx1 * Lz1
                                   + buff[CONST_NDIM_2*icmp + 1] * Lx2 * Lz1
                                   + buff[CONST_NDIM_2*icmp + 2] * Lx1 * Lz2
                                   + buff[CONST_NDIM_2*icmp + 3] * Lx2 * Lz2;
    }
  }
  CUDACHECK(cudaFree(buff_d));
  CUDACHECK(cudaFree(indx1d_d));

  return 0;
}

int
io_line_keep(ioline_t *ioline, float * w_end_d,
             float * buff, int it, int ncmp, int siz_icmp)
{
  int size = sizeof(float)*ncmp;
  float *buff_d = (float *) cuda_malloc(size);
  dim3 block(32);
  dim3 grid;
  grid.x = (ncmp+block.x-1)/block.x;

  for (int n=0; n < ioline->num_of_lines; n++)
  {
    int   *this_line_iptr   = ioline->recv_iptr[n];
    float *this_line_seismo = ioline->recv_seismo[n];
  
    for (int ir=0; ir < ioline->line_nr[n]; ir++)
    {
      int iptr = this_line_iptr[ir];
      float *this_seismo = this_line_seismo + ir * ioline->max_nt * ncmp;
      io_recv_line_pack_buff<<<grid, block>>>(w_end_d, buff_d, ncmp, siz_icmp, iptr);
      CUDACHECK(cudaMemcpy(buff,buff_d,size,cudaMemcpyDeviceToHost));
      for (int icmp=0; icmp < ncmp; icmp++)
      {
        int iptr_seismo = icmp * ioline->max_nt + it;
        this_seismo[iptr_seismo] = buff[icmp];
      }
    }
  }
  CUDACHECK(cudaFree(buff_d));
  return 0;
}

__global__ void
io_recv_line_interp_pack_buff(float *var, float *buff_d, int ncmp, size_t siz_icmp, size_t *indx1d_d)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  //indx1d_d size is CONST_NDIM_2 = 4
  if(ix < ncmp)
  {
   buff_d[4*ix+0] = var[ix*siz_icmp + indx1d_d[0] ];
   buff_d[4*ix+1] = var[ix*siz_icmp + indx1d_d[1] ];
   buff_d[4*ix+2] = var[ix*siz_icmp + indx1d_d[2] ];
   buff_d[4*ix+3] = var[ix*siz_icmp + indx1d_d[3] ];
  }
}

__global__ void
io_recv_line_pack_buff(float *var, float *buff_d, int ncmp, size_t siz_icmp, int iptr)
{
  size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
  if(ix < ncmp)
  {
   buff_d[ix] = var[ix*siz_icmp + iptr];
  }
}



int
io_recv_output_sac(iorecv_t *iorecv,
                   float dt,
                   int num_of_vars,
                   char **cmp_name,
                   char *evtnm,
                   char *output_dir,
                   char *err_message)
{
  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;
  char ou_file[CONST_MAX_STRLEN];

  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;

    //fprintf(stdout,"=== Debug: num_of_vars=%d\n",num_of_vars);fflush(stdout);
    for (int icmp=0; icmp < num_of_vars; icmp++)
    {
      //fprintf(stdout,"=== Debug: icmp=%d\n",icmp);fflush(stdout);

      float *this_trace = this_recv->seismo + icmp * iorecv->max_nt;

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm,
                      this_recv->name, cmp_name[icmp]);

      //fprintf(stdout,"=== Debug: icmp=%d,ou_file=%s\n",icmp,ou_file);fflush(stdout);

      sacExport1C1R(ou_file,
            this_trace,
            evt_x, evt_y, evt_z, evt_d,
            this_recv->x, 0.0, this_recv->z,
            dt, dt, iorecv->max_nt, err_message);
    }
  }

  return 0;
}

int io_line_output_sac(ioline_t *ioline,
      float dt, char **cmp_name, char *evtnm, char *output_dir)
{
  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;
  char ou_file[CONST_MAX_STRLEN];
  char err_message[CONST_MAX_STRLEN];
  
  for (int n=0; n < ioline->num_of_lines; n++)
  {
    int   *this_line_iptr   = ioline->recv_iptr[n];
    float *this_line_seismo = ioline->recv_seismo[n];

    for (int ir=0; ir < ioline->line_nr[n]; ir++)
    {
      float *this_seismo = this_line_seismo + ir * ioline->max_nt * ioline->ncmp;

      for (int icmp=0; icmp < ioline->ncmp; icmp++)
      {
        float *this_trace = this_seismo + icmp * ioline->max_nt;

        sprintf(ou_file,"%s/%s.%s.no%d.%s.sac", output_dir,evtnm,
                  ioline->line_name[n],ioline->recv_seq[n][ir],
                  cmp_name[icmp]);

        sacExport1C1R(ou_file,
              this_trace,
              evt_x, evt_y, evt_z, evt_d,
              ioline->recv_x[n][ir],
              0.0,
              ioline->recv_z[n][ir],
              dt, dt, ioline->max_nt, err_message);
      } // icmp
    } // ir
  } // line

  return 0;
}

// calculate and output strain cmp for elastic medium
//   do not find a better file to hold this func
//   temporarily put here

int
io_recv_output_sac_el_iso_strain(iorecv_t *iorecv,
                     float * lam3d,
                     float * mu3d,
                     float dt,
                     char *evtnm,
                     char *output_dir,
                     char *err_message)
{
  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;
  char ou_file[CONST_MAX_STRLEN];

  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;
    int iptr = this_recv->indx1d[0];

    float lam = lam3d[iptr];
    float mu  =  mu3d[iptr];

    // cmp seq hard-coded, need to revise in the future
    float *Txx = this_recv->seismo + 2 * iorecv->max_nt;
    float *Tzz = this_recv->seismo + 3 * iorecv->max_nt;
    float *Txz = this_recv->seismo + 4 * iorecv->max_nt;

    float E1 = (lam + mu) / (mu * ( 3.0 * lam + 2.0 * mu));
    float E2 = - lam / ( 2.0 * mu * (3.0 * lam + 2.0 * mu));
    float E3 = 1.0 / mu;

    // conver to strain per time step
    for (int it = 0; it < iorecv->max_nt; it++)
    {
      float E0 = E2 * (Txx[it] + Tzz[it]);

      Txx[it] = E0 - (E2 - E1) * Txx[it];
      Tzz[it] = E0 - (E2 - E1) * Tzz[it];
      Txz[it] = 0.5 * E3 * Txz[it];
    }

    // output to sca file
    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exx");
    sacExport1C1R(ou_file,Txx,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, 0.0, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Ezz");
    sacExport1C1R(ou_file,Tzz,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, 0.0, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exz");
    sacExport1C1R(ou_file,Txz,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, 0.0, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

  } // loop ir

  return 0;
}

int iosnap_print(iosnap_t *iosnap)
{    
  fprintf(stdout, "--> num_of_snap = %d\n", iosnap->num_of_snap);
  fprintf(stdout, "#   i0 k0 ni nk di dk it0 dit vel stress strain gi1 gk1\n");
  for (int n=0; n < iosnap->num_of_snap; n++)
  {
    fprintf(stdout, " %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
              n,
              iosnap->i1[n], iosnap->k1[n],
              iosnap->ni[n], iosnap->nk[n],
              iosnap->di[n], iosnap->dk[n],
              iosnap->it1[n], iosnap->dit[n], 
              iosnap->out_vel[n],
              iosnap->out_stress[n],
              iosnap->out_strain[n],
              iosnap->i1_to_glob[n],
              iosnap->k1_to_glob[n]);
  }

  return 0;
}

int iorecv_print(iorecv_t *iorecv)
{    
  //fprintf(stdout, "\n");
  //fprintf(stdout, "--> station information.\n");
  //fprintf(stdout, " number_of_station  = %4d\n", blk->number_of_station);
  //fprintf(stdout, " seismo_format_sac  = %4d\n", blk->seismo_format_sac );
  //fprintf(stdout, " seismo_format_segy = %4d\n", blk->seismo_format_segy);
  //fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  //fprintf(stdout, "\n");

  //if(blk->number_of_station > 0)
  //{
  //    //fprintf(stdout, " station_indx:\n");
  //    fprintf(stdout, " stations             x           z           i           k:\n");
  //}

  //for(n=0; n<blk->number_of_station; n++)
  //{
  //    indx = 2*n;
  //    fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //            blk->station_coord[indx], blk->station_coord[indx+1],
  //            blk->station_indx [indx], blk->station_indx [indx+1]);
  //}
  //fprintf(stdout, "\n");

  return 0;
}

/*
 * get next non-comment line
 */

int
io_get_nextline(FILE *fp, char *str, int length)
{
  int ierr = 0;

  do
  {
    if (fgets(str, length, fp) == NULL)
    {
      ierr = 1;
      return ierr;
    }
  } while (str[0] == '#' || str[0] == '\n');

  // remove newline char
  int len = strlen(str);

  if (len > 0 && str[len-1] == '\n') {
    str[len-1] = '\0';
  }

  // for debug:
  //fprintf(stdout," --return: %s\n", str);

  return ierr;
}
