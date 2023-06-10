/*
 * source term related processing
 */

// todo:

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "src_t.h"
#include "io_funcs.h"

/*
 * src_t alloc
 */

int
src_init(src_t *src, int force_actived, int moment_actived,
         int num_of_src, int max_nt, int max_stage, int max_ext)
{
  // set default value
  src->total_number = num_of_src;
  src->max_nt    = max_nt;
  src->max_stage = max_stage;
  src->max_ext   = max_ext;

  src->force_actived   = force_actived;
  src->moment_actived   = moment_actived;

  // allocate var
  src->si = (int *)malloc(num_of_src*sizeof(int));
  src->sk = (int *)malloc(num_of_src*sizeof(int));
  src->it_begin = (int *)malloc(num_of_src*sizeof(int));
  src->it_end   = (int *)malloc(num_of_src*sizeof(int));
  src->ext_num  = (int *)malloc(num_of_src*sizeof(int));
  src->ext_indx = (int *)malloc(num_of_src*max_ext * sizeof(int  ));
  src->ext_coef = (float *)malloc(num_of_src*max_ext * sizeof(float));

  src->Fx = NULL;
  src->Fz = NULL;

  if (force_actived == 1) {
    src->Fx = (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Fz = (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    for (int iptr=0; iptr < max_stage * max_nt * num_of_src; iptr++) {
      src->Fx[iptr] = 0.0;
      src->Fz[iptr] = 0.0;
    }
  }

  src->Mxx = NULL;
  src->Mzz = NULL;
  src->Mxz = NULL;

  if (moment_actived == 1) {
    src->Mxx= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mzz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mxz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    for (int iptr=0; iptr < max_stage * max_nt * num_of_src; iptr++) {
      src->Mxx[iptr] = 0.0;
      src->Mzz[iptr] = 0.0;
      src->Mxz[iptr] = 0.0;
    }
  }

  return 0;
}

int
src_set_time(src_t *src, int it, int istage)
{
  src->it     = it;
  src->istage = istage;

  return 0;
}

/*
 * read .src file and convert into internal structure
 */

int
src_read_locate_file(gd_t     *gd,
                     src_t    *src,
                     char     *in_src_file,
                     float     t0,
                     float     dt,
                     int       max_stage,
                     float    *rk_stage_time,
                     int       npoint_half_ext,
                     int       verbose)
{
  int ierr = 0;

  // get grid info from gd
  int   ni1 = gd->ni1;
  int   ni2 = gd->ni2;
  int   nk1 = gd->nk1;
  int   nk2 = gd->nk2;
  int   nx  = gd->nx ;
  int   nz  = gd->nz ;
  int   npoint_ghosts = gd->npoint_ghosts;
  size_t siz_iz= gd->siz_iz;

  // get total elem of exted src region for a single point
  int len_ext = 2*npoint_half_ext+1;
  int max_ext = len_ext * len_ext;

  // local
  FILE *fp =NULL;
  char str[500];

  // numbe of source, could be force and/or moment
  int in_num_source;
  // input location is grid index (0) or coordinate (1)
  int is_coord;
  // the 2rd coord z is coord (0) or depth (1)
  int is_depth;
  // stf is specified by wavelet name (0) or values (1)
  int in_stf_given;
  float in_stf_length=0.0;
  float in_stf_dt=0.0;
  int   in_stf_nt=1;
  // which cmp is used, 1: force; 2: moment, 3: force + moment
  int in_cmp_type;
  // moment is given by tensor (0) or angle + mu D A (1)
  int in_mechanism_type;

  // open in_src_file
  if ((fp = fopen(in_src_file, "r"))==NULL) {
    fprintf(stderr,"ERROR: fail to open in_src_file=%s", in_src_file);
    fflush(stderr); exit(1);
  }

  // event name
  if (!io_get_nextline(fp, str,500)) {
    sprintf(src->evtnm,"%s",str);
  }
  // number of source
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d",&in_num_source);
  }
  if (in_num_source <= 0) {
    fprintf(stderr,"ERROR: in_num_source=%d <=0\n", in_num_source);
    fflush(stderr); exit(1);
  }

  // source time function is given by wavelet name or values
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d",&in_stf_given);
    if (in_stf_given == 0) { // by name
      sscanf(str,"%d %f",&in_stf_given, &in_stf_length);
    } else if (in_stf_given == 1) { // by value
      sscanf(str,"%d %f %d",&in_stf_given, &in_stf_dt, &in_stf_nt);
    } else {
      fprintf(stderr, "ERROR: in_stf_given=%d invalid (either 0 or 1)\n", in_stf_given);
      fflush(stderr); exit(1);
    }
  }

  // force and/or moment, and moment by tensor or angle + muDa
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d %d",&in_cmp_type, &in_mechanism_type);
  }

  // meaning of location and the 2rd input if location is given by coord
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d %d",&is_coord, &is_depth);
  }

  //
  // loop each source to get locate and index
  //

  int num_of_src_here = 0;

  float sx, sz;
  int   si, sk;
  float sx_inc, sz_inc;

  float **all_inc   = fdlib_mem_calloc_2l_float(in_num_source,CONST_NDIM,0.0,"all_inc");
  int   **all_index = fdlib_mem_calloc_2l_int(in_num_source,CONST_NDIM,0,"all_index");
  int   *all_in_thread = fdlib_mem_calloc_1d_int(in_num_source,0,"all_in_thread");

  // read coords and determine if in this thread
  for (int is=0; is<in_num_source; is++)
  {
    // read in and get src global index
    if (!io_get_nextline(fp, str,500))
    {
      // read in as float value
      sscanf(str,"%f %f", &sx, &sz);

      if (is_coord == 1) // physical coord
      {
        // convert to global index
        //    todo: check if out of computational region
        if (gd->type == GD_TYPE_CURV)
        {
          // if sz is depth, convert to axis when it is in this thread
          if (is_depth == 1) {
            gd_curv_depth_to_axis(gd,sx,&sz);
          }
          gd_curv_coord_to_local_indx(gd,sx,sz,
                                      &si,&sk,&sx_inc,&sz_inc);
        }
        else if (gd->type == GD_TYPE_CART)
        {
          // if sz is depth, convert to axis
          if (is_depth == 1) {
            sz = gd->z1d[gd->nk2] - sz;
          }
          gd_cart_coord_to_local_indx(gd,sx,sz,
                                      &si,&sk,&sx_inc,&sz_inc);
        }
        // keep index to avoid duplicat run
        all_index[is][0] = si;
        all_index[is][1] = sk;
        all_inc  [is][0] = sx_inc;
        all_inc  [is][1] = sz_inc;

        //-- to notice user the progress using screen output for large input
        if ((is % 1000 ==0) && verbose>99) {
          fprintf(stdout,"-- loc %d-th src index, finish %2.0f%%\n",
                      is, (float)(is+1)/in_num_source*100.0);
          fflush(stdout);
        }
      }
      else // computational coordinate or grid index
      {
        // add ghosts, to local point
        sx = sx + gd->fdx_nghosts;
        sz = sz + gd->fdz_nghosts;
        // if sz is relative to surface, convert to normal index
        if (is_depth == 1) {
          sz = gd->nk2 - sz + gd->fdz_nghosts;
        }

        // nearest integer index
        si = (int) (sx + 0.5);
        sk = (int) (sz + 0.5);
        // relative shift
        sx_inc = sx - si;
        sz_inc = sz - sk;

        all_index[is][0] = si;
        all_index[is][1] = sk;
        all_inc  [is][0] = sx_inc;
        all_inc  [is][1] = sz_inc;
      }
    }

    // check if in this thread using index
    if (gd_lindx_is_inner(si,sk,gd)==1)
    {
      num_of_src_here += 1;
      all_in_thread[is] = 1;
    }
  } // is loop

  // print for QC
  if (verbose > 99) {
    fprintf(stdout,"src located results:\n");
    for (int is=0; is < in_num_source; is++)
    {
      fprintf(stdout,"-- %d: indx=(%d,%d), inc=(%f,%f)\n",
                    is, all_index[is][0],all_index[is][1],
                        all_inc[is][0],all_inc[is][1]);
    }
    fflush(stdout);
  }

  if (verbose > 1) {
    for (int is=0; is < in_num_source; is++)
    {
      if(all_in_thread[is] == 0)
      {
        fprintf(stdout,"#########         ########\n");
        fprintf(stdout,"######### Warning ########\n");
        fprintf(stdout,"#########         ########\n");
        fprintf(stdout,"source number %d physical coordinates are outside calculation area !\n",is);
      }
    }
    fflush(stdout);
  }

  //
  // alloc src_t struct for this thread
  //

  // check if force and moment used
  int force_actived  = 0;
  int moment_actived = 0;
  if (num_of_src_here > 0)
  {
    if (in_cmp_type == 1 || in_cmp_type == 3) {
      force_actived = 1;
    }

    if (in_cmp_type == 2 || in_cmp_type == 3) {
      moment_actived = 1;
    } 
  }

  // get number of sample for src_t
  int max_nt = 0;
  if (in_stf_given == 0) { // by name
    max_nt = (int) (in_stf_length / dt + 0.5);
  } else { // by value
    max_nt = (int) (((in_stf_nt-1)*in_stf_dt / dt)+ 0.5) + 1; 
  }

  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  //
  // loop all source and only keep those in this thread
  //

  float wavelet_tstart;
  char  wavelet_name[50]; // assume max size of name is <= 50
  float wavelet_coefs[10]; // assume max number of coef <= 10
  int  it_begin;
  int  it_end;

  float fx,fz;
  float mxx,mzz,mxz;
  float *f1 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"f1");
  float *f3 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"f3");
  float *m11 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m11");
  float *m33 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m33");
  float *m13 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m13");
  float *t_in = (float *)malloc(in_stf_nt*sizeof(float));

  int is_local = 0;
  for (int is=0; is<in_num_source; is++)
  {
    // read stf and cmp of each source
    if (in_stf_given == 0) // wavelet name
    {
      // read in stf
      if (!io_get_nextline(fp, str,500))
      {
        if (all_in_thread[is] == 1) // in in this thread
        {
          // read up to 10 coefs, may be less than 10
          sscanf(str,"%f %s %f %f %f %f %f %f %f %f %f %f",
                  &wavelet_tstart, wavelet_name, wavelet_coefs+0,
                  wavelet_coefs+1, wavelet_coefs+2, wavelet_coefs+3,
                  wavelet_coefs+4, wavelet_coefs+5, wavelet_coefs+6,
                  wavelet_coefs+7, wavelet_coefs+8, wavelet_coefs+9);
        }
      }
      // read in cmp
      if (!io_get_nextline(fp, str,500))
      {
        if (all_in_thread[is] == 1) // in in this thread
        {
          if (in_cmp_type == 1) { // force
            sscanf(str,"%f %f ",&fx,&fz);
          } else if (in_cmp_type == 2) { // moment
            sscanf(str,"%f %f %f ",&mxx,&mzz,&mxz);
          } else { // force + moment
            sscanf(str,"%f %f %f %f %f ", &fx,&fz,&mxx,&mzz,&mxz);
          }
        }
      }
    }
    else // by values
    {
      // read t0
      if (!io_get_nextline(fp, str,500)) {
        sscanf(str,"%f",&wavelet_tstart);
      }

      // read cmp in number of in_stf_nt no matter in_thread or not
      for (int it=0; it<in_stf_nt; it++)
      {
        if (!io_get_nextline(fp, str,500))
        {
          if (all_in_thread[is] == 1) // in in this thread
          {
            if (in_cmp_type == 1) { // force
              sscanf(str,"%f %f",f1+it,f3+it);
            } else if (in_cmp_type == 2) { // moment
              sscanf(str,"%f %f %f",m11+it,m33+it,m13+it);
            } else { // force + moment
              sscanf(str,"%f %f %f %f %f",
                  f1+it,f3+it,m11+it,m33+it,m13+it);
            }
          } // in this thread
        } // get next line
      } // it
    } // read in stf for one is

    // push into src_t if in this thread
    if (all_in_thread[is] == 1)
    {
      si = all_index[is][0];
      sk = all_index[is][1];
  
      // keep into src_t
      src->si[is_local] = si;
      src->sk[is_local] = sk;

      // for extended points and coefs
      sx_inc = all_inc[is][0];
      sz_inc = all_inc[is][1];
      float wid_gauss = npoint_half_ext / 2.0;
      float *this_ext_coef = src->ext_coef + is_local * max_ext;

      src_cal_norm_delt2d(this_ext_coef, sx_inc, sz_inc,
                          wid_gauss, wid_gauss, npoint_half_ext);

      size_t iptr_ext = 0;
      for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
      {
        for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
        {
          if (gd_lindx_is_inner(i,k,gd)==1)
          {
            // Note index need match coef
            int iptr_grid = i + k * siz_iz;
            int iptr_coef = (i-(si-npoint_half_ext))
                            + len_ext * (k-(sk-npoint_half_ext)); 
            src->ext_indx[iptr_ext + is_local * max_ext] = iptr_grid;
            src->ext_coef[iptr_ext + is_local * max_ext] = this_ext_coef[iptr_coef];
            iptr_ext++;
          }
        }
      }
      // only count index inside phys region for this thread
      src->ext_num[is_local] = iptr_ext;

      //
      // wavelet
      //

      // time step, considering t0
      it_begin = (int) ( (wavelet_tstart - t0) / dt);
      it_end   = it_begin + max_nt - 1;

      src->it_begin[is_local] = it_begin;
      src->it_end  [is_local] = it_end  ;

      // setp input t vector for interp 
      for(int it=0; it<in_stf_nt; it++)
      {
        t_in[it] = wavelet_tstart + it*in_stf_dt;
      }

      for (int it=it_begin; it<=it_end; it++)
      {
        int it_to_it1 = (it - it_begin);
        int iptr_it = is_local * max_nt * max_stage + it_to_it1 * max_stage;
        // need to explain purpose
        float t_shift = wavelet_tstart - (it_begin * dt + t0);

        for (int istage=0; istage<max_stage; istage++)
        {
          int iptr = iptr_it + istage;

          // cal stf for given wavelet name
          if (in_stf_given==0)
          {
            // time relative to start time of this source, considering diff from int conversion
            float t = it_to_it1 * dt + rk_stage_time[istage] * dt - t_shift;

            float stf_val = src_cal_wavelet(t,wavelet_name,wavelet_coefs);

            if (force_actived==1) {
              src->Fx[iptr]  = stf_val * fx;
              src->Fz[iptr]  = stf_val * fz;
            }
            if (moment_actived==1) {
              src->Mxx[iptr] = stf_val * mxx;
              src->Mzz[iptr] = stf_val * mzz;
              src->Mxz[iptr] = stf_val * mxz;
            }
          }
          // interp for input values
          else
          {
            // time relative to start time of this source, considering diff from int conversion
            float t = it * dt + rk_stage_time[istage] * dt - t_shift;

            // interp1d order
            int order = 3;     

            if (force_actived==1)
            {
              fx = LagInterp_Piecewise_1d(t_in, f1, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              fz = LagInterp_Piecewise_1d(t_in, f3, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);

              src->Fx[iptr]  = fx;
              src->Fz[iptr]  = fz;
            }

            if (moment_actived==1)
            {
              mxx = LagInterp_Piecewise_1d(t_in, m11, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              mzz = LagInterp_Piecewise_1d(t_in, m33, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              mxz = LagInterp_Piecewise_1d(t_in, m13, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              src->Mxx[iptr] = mxx;
              src->Mzz[iptr] = mzz;
              src->Mxz[iptr] = mxz;
            }
          }
        } // istage
      } // it
      
      // local is increase
      is_local += 1;

    } // if in_thread
  } // is

  // close file and free local pointer
  fclose(fp); 

  free(f1);
  free(f3);
  free(m11);
  free(m33);
  free(m13);
  free(t_in);
  free(all_in_thread);

  fdlib_mem_free_2l_float(all_inc, in_num_source, "free_all_inc");
  fdlib_mem_free_2l_int  (all_index, in_num_source, "free_all_index");

  return ierr;
}

/*
 * 2d spatial smoothing
 */

int
src_cal_norm_delt2d(float *delt, float x0, float z0,
                    float rx0, float rz0, int LenDelt)
{
  float SUM = 0.0 ;

  int iptr = 0;
  for(int k=-LenDelt; k<=LenDelt; k++) {
      for(int i=-LenDelt; i<=LenDelt; i++) {
        float D1 = fun_gauss(i-x0, rx0 ,0.0);           
        float D3 = fun_gauss(k-z0, rz0 ,0.0);          
        delt[iptr] = D1 * D3;
        SUM += delt[iptr];
        iptr++;
      }
  }

  if( SUM < 1e-20 )
  {
    fprintf(stderr, "cal_norm_delt is zero\n");
    exit(1);
  }

  int siz_1d = 2 * LenDelt + 1;
  for (int iptr=0; iptr< siz_1d*siz_1d; iptr++) {
    delt[iptr] /= SUM;
  }

  return 0;
} 

/*
 * get stf value at a given t
 */

float
src_cal_wavelet(float t, char *wavelet_name, float *wavelet_coefs)
{
  float stf_val;

  if (strcmp(wavelet_name, "ricker")==0) {
    stf_val = fun_ricker(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "gaussian")==0) {
    stf_val = fun_gauss(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "ricker_deriv")==0) {
    stf_val = fun_ricker_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "gaussian_deriv")==0) {
    stf_val = fun_gauss_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else{
    fprintf(stderr,"wavelet_name=%s\n", wavelet_name); 
    fprintf(stderr,"   not implemented yet\n"); 
    fflush(stderr); exit(1);
  }

  return stf_val;
}

/*
 * wavelet functions
 */

// ricker and it deriv.
float 
fun_ricker(float t, float fc, float t0)
{
  float u = (t-t0)*2.0*PI*fc;
  float v = (1-u*u/2.0)*exp(-u*u/4.0);

  return v;
}

float 
fun_ricker_deriv(float t, float fc, float t0)
{
  float u = (t-t0)*2.0*PI*fc;
  float v = u*(-3+1/2*u*u)*exp(-u*u/4)*PI*fc;

  return v;
}
//gauss and it deriv
float
fun_gauss(float t, float a, float t0)
{
  float f;
  f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrtf(PI)*a);
  return f;
}

float
fun_gauss_deriv(float t, float a, float t0)
{
  float f;
  f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrtf(PI)*a)*(-2*(t-t0)/(a*a));
  return f;
}


