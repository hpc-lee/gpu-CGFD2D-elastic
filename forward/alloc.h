#ifndef ALLOC_H
#define ALLOC_H

#include "constants.h"
#include "gd_t.h"
#include "fd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"

int
init_gd_device(gd_t *gd, gd_t *gd_d);

int
init_md_device(md_t *md, md_t *md_d);

int
init_fd_device(fd_t *fd, fd_wav_t *fd_wav_d);

int
init_metric_device(gdcurv_metric_t *metric, gdcurv_metric_t *metric_d);

int
init_src_device(src_t *src, src_t *src_d);

int 
init_bdry_device(gd_t *gd, bdry_t *bdry, bdry_t *bdry_d);

int 
init_wave_device(wav_t *wav, wav_t *wav_d);

int 
dealloc_md_device(md_t md_d);

int 
dealloc_fd_device(fd_wav_t fd_wav_d);

int
dealloc_metric_device(gdcurv_metric_t metric_d);

int 
dealloc_src_device(src_t src_d);

int 
dealloc_bdry_device(bdry_t bdry_d);

int 
dealloc_wave_device(wav_t wav_d);

#endif
