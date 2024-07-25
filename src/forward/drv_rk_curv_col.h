#ifndef DRV_RK_CURV_COL_H
#define DRV_RK_CURV_COL_H

#include "fd_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

int
drv_rk_curv_col_allstep(
  fd_t        *fd,
  gd_t        *gd,
  gdcurv_metric_t *metric,
  md_t      *md,
  src_t      *src,
  bdry_t *bdry,
  wav_t  *wav,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose);

#endif
