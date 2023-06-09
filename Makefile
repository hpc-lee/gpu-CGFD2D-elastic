################################################################################
#  Makefile for CGFD2D-wave package
#
#  Author: ZHANG Wei <zhangwei@sustech.edu.cn>
#  Copyright (C) Wei ZHANG, 2021. All Rights Reserved.
#
################################################################################

# known issues:
#  .h dependency is not included, use make cleanall

#-------------------------------------------------------------------------------
# compiler
#-------------------------------------------------------------------------------

CXX    := $(GNUHOME)/bin/g++
GC     := $(CUDAHOME)/bin/nvcc 

#- debug
#CFLAGS   := -g -std=c99 $(CFLAGS)
#CPPFLAGS := -g -std=c++11 $(CPPFLAGS)
#- O3
CPPFLAGS := -O3 -std=c++11 $(CPPFLAGS)

CFLAGS_CUDA   := -O3 -arch=sm_80 -std=c++11 -w -rdc=true
CFLAGS_CUDA += -I$(CUDAHOME)/include
CFLAGS_CUDA += -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/
#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf
LDFLAGS += -lm -arch=sm_80

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

main_curv_col_2d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o interp.o alloc.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry2d.o \
		media_read_file.o \
		gd_t.o md_t.o wav_t.o \
		bdry_t.o src_t.o io_funcs.o \
		blk_t.o cuda_common.o\
		drv_rk_curv_col.o \
		sv_curv_col_el_gpu.o \
		sv_curv_col_ac_iso_gpu.o \
		sv_curv_col_el_iso_gpu.o \
		sv_curv_col_el_vti_gpu.o \
		sv_curv_col_el_aniso_gpu.o \
		main_curv_col_2d.o
	$(GC) -o $@ $^ $(LDFLAGS)

media_geometry2d.o: media/media_geometry2d.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_utility.o: media/media_utility.cpp 
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_layer2model.o: media/media_layer2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_grid2model.o: media/media_grid2model.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
media_read_file.o: media/media_read_file.cpp
	${CXX} -c -o $@ $(CPPFLAGS) $<
cJSON.o: lib/cJSON.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sacLib.o: lib/sacLib.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
fdlib_mem.o: lib/fdlib_mem.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
fdlib_math.o: lib/fdlib_math.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
fd_t.o: forward/fd_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
par_t.o: forward/par_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
alloc.o: forward/alloc.cu
	$(GC) -c -o $@ $(CFLAGS_CUDA) $<
cuda_common.o: forward/cuda_common.cu
	$(GC) -c -o $@ $(CFLAGS_CUDA) $<
interp.o: forward/interp.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
gd_t.o: forward/gd_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
md_t.o: forward/md_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
wav_t.o: forward/wav_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
bdry_t.o: forward/bdry_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
src_t.o: forward/src_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
io_funcs.o: forward/io_funcs.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
blk_t.o: forward/blk_t.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
drv_rk_curv_col.o:          forward/drv_rk_curv_col.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sv_curv_col_el_gpu.o:          forward/sv_curv_col_el_gpu.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sv_curv_col_el_iso_gpu.o:   forward/sv_curv_col_el_iso_gpu.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sv_curv_col_el_vti_gpu.o:   forward/sv_curv_col_el_vti_gpu.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sv_curv_col_el_aniso_gpu.o: forward/sv_curv_col_el_aniso_gpu.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<
sv_curv_col_ac_iso_gpu.o:   forward/sv_curv_col_ac_iso_gpu.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<

main_curv_col_2d.o: forward/main_curv_col_2d.cu
	${GC} -c -o $@ $(CFLAGS_CUDA) $<

cleanexe:
	rm -f main_curv_col_2d

cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
