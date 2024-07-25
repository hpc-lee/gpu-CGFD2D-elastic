#ifndef IO_FUNCS_H
#define IO_FUNCS_H

#include "constants.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"

/*************************************************
 * structure
 *************************************************/

// for stations output

// single station
typedef struct
{
  float x;
  float z;
  float di;
  float dk;
  int   i;
  int   k;
  size_t   indx1d[CONST_NDIM_2];
  float *seismo;
  char  name[CONST_MAX_STRLEN];
} iorecv_one_t;

typedef struct
{
  int                 total_number;
  int                 max_nt;
  int                 ncmp;
  iorecv_one_t *recvone;
} iorecv_t;

// line output
typedef struct
{
  int     num_of_lines; 
  int     max_nt;
  int     ncmp;

  int    *line_nr; // number of receivers, for name from input file
  int    *line_seq; // line number, for name from input file
  //int    **recv_ir;
  //int    **recv_jr;
  //int    **recv_kr;
  int    **recv_seq; // recv seq in this line
  int    **recv_iptr;
  float  **recv_x; // for sac output
  float  **recv_z; // for sac output
  float  **recv_seismo;
  char   **line_name;
} ioline_t;

// snapshot output
typedef struct
{
  // for esti size of working space var
  size_t siz_max_wrk;

  int num_of_snap;

  int *i1;
  int *k1;
  int *ni;
  int *nk;
  int *di;
  int *dk;
  int *it1;
  int *dit;
  int *out_vel;
  int *out_stress;
  int *out_strain;

  int *i1_to_glob;
  int *k1_to_glob;

  char **fname;
} iosnap_t;

// for nc output

typedef struct
{
  int num_of_snap;
  int *ncid;
  int *timeid;
  int *varid_V;  // [num_of_snap*CONST_NDIM];
  int *varid_T;  // [num_of_snap*CONST_NDIM_2];
  int *varid_E;  // [num_of_snap*CONST_NDIM_2];
  int *cur_it ;  // [num_of_snap];
}
iosnap_nc_t;

/*************************************************
 * function prototype
 *************************************************/

int
io_recv_read_locate(gd_t *gd,
                    iorecv_t  *iorecv,
                    int       nt_total,
                    int       num_of_vars,
                    char *in_filenm);

int
io_line_locate(gd_t *gd,
               ioline_t *ioline,
               int    num_of_vars,
               int    nt_total,
               int    number_of_receiver_line,
               int   *receiver_line_index_start,
               int   *receiver_line_index_incre,
               int   *receiver_line_count,
               char **receiver_line_name);

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
                    char *output_dir);

int
io_snap_nc_create(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc);

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
               float time);

int
io_snap_nc_create_ac(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc);

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
               float time);

__global__ void
io_snap_pack_buff(float * var,
                  size_t siz_iz,
                  int starti,
                  int counti,
                  int increi,
                  int startk,
                  int countk,
                  int increk,
                  float * buff_d);

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
                               int increk);
int
io_snap_nc_close(iosnap_nc_t *iosnap_nc);

int
io_recv_keep(iorecv_t *iorecv, float * w_end_d,
             float * buff, int it, int ncmp, int siz_icmp);

int
io_line_keep(ioline_t *ioline, float * w_end_d,
             float * buff, int it, int ncmp, int siz_icmp);

__global__ void
io_recv_line_interp_pack_buff(float *var, float *buff_d, int ncmp, size_t siz_icmp, size_t *indx1d_d);

__global__ void
io_recv_line_pack_buff(float *var, float *buff_d, int ncmp, size_t siz_icmp, int iptr);

int
io_recv_output_sac(iorecv_t *iorecv,
                   float dt,
                   int num_of_vars,
                   char **cmp_name,
                   char *evtnm,
                   char *output_dir,
                   char *err_message);

int
io_line_output_sac(ioline_t *ioline,
      float dt, char **cmp_name, char *evtnm, char *output_dir);

int
io_recv_output_sac_el_iso_strain(iorecv_t *iorecv,
                     float * lam3d,
                     float * mu3d,
                     float dt,
                     char *evtnm,
                     char *output_dir,
                     char *err_message);

int
iosnap_print(iosnap_t *iosnap);

int
iorecv_print(iorecv_t *iorecv);

int
io_get_nextline(FILE *fp, char *str, int length);

#endif
