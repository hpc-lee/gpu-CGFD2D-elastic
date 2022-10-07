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

CC     := /data3/lihl/software/gcc-10.3.0-compile/bin/gcc
CXX    := /data3/lihl/software/gcc-10.3.0-compile/bin/g++
NETCDF := /data3/lihl/software/disable-netcdf-4.4.1

#-- 
CFLAGS := -I$(NETCDF)/include -I./lib/ -I./forward/ -I./media/  $(CFLAGS)

#- debug
#CFLAGS   := -g -std=c99 $(CFLAGS)
#CPPFLAGS := -g -std=c++11 $(CPPFLAGS)
#- O3
CFLAGS   := -O3 -std=c99 $(CFLAGS)
CPPFLAGS := -std=c++11 $(CPPFLAGS)

#- static
#LDFLAGS := $(NETCDF)/lib/libnetcdf.a -lm -static $(LDFLAGS)
#LDFLAGS := -lm  $(LDFLAGS) $(NETCDF)/lib/libnetcdf.a
#- dynamic
LDFLAGS := -L$(NETCDF)/lib -lnetcdf -lm $(LDFLAGS)

#-------------------------------------------------------------------------------
# target
#-------------------------------------------------------------------------------

# special vars:
# 	$@ The file name of the target
# 	$< The names of the first prerequisite
#   $^ The names of all the prerequisites 

main_curv_col_2d: \
		cJSON.o sacLib.o fdlib_mem.o fdlib_math.o  \
		fd_t.o par_t.o \
		media_utility.o \
		media_layer2model.o \
		media_grid2model.o \
		media_geometry2d.o \
		media_read_file.o \
		gd_t.o md_t.o wav_t.o \
		bdry_t.o src_t.o io_funcs.o \
		blk_t.o interp.o\
		drv_rk_curv_col.o \
		sv_curv_col_el.o \
		sv_curv_col_ac_iso.o \
		sv_curv_col_el_iso.o \
		sv_curv_col_el_vti.o \
		sv_curv_col_el_aniso.o \
		main_curv_col_2d.o
	$(CXX) -o $@ $^ $(LDFLAGS)

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
cJSON.o: lib/cJSON.c
	${CC} -c -o $@ $(CFLAGS) $<
sacLib.o: lib/sacLib.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_mem.o: lib/fdlib_mem.c
	${CC} -c -o $@ $(CFLAGS) $<
fdlib_math.o: lib/fdlib_math.c
	${CC} -c -o $@ $(CFLAGS) $<
fd_t.o: forward/fd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
par_t.o: forward/par_t.c
	${CC} -c -o $@ $(CFLAGS) $<
interp.o: forward/interp.c
	${CC} -c -o $@ $(CFLAGS) $<
gd_t.o: forward/gd_t.c
	${CC} -c -o $@ $(CFLAGS) $<
md_t.o: forward/md_t.c
	${CC} -c -o $@ $(CFLAGS) $<
wav_t.o: forward/wav_t.c
	${CC} -c -o $@ $(CFLAGS) $<
bdry_t.o: forward/bdry_t.c
	${CC} -c -o $@ $(CFLAGS) $<
src_t.o: forward/src_t.c
	${CC} -c -o $@ $(CFLAGS) $<
io_funcs.o: forward/io_funcs.c
	${CC} -c -o $@ $(CFLAGS) $<
blk_t.o: forward/blk_t.c
	${CC} -c -o $@ $(CFLAGS) $<
drv_rk_curv_col.o:          forward/drv_rk_curv_col.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_curv_col_el.o:          forward/sv_curv_col_el.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_curv_col_el_iso.o:   forward/sv_curv_col_el_iso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_curv_col_el_vti.o:   forward/sv_curv_col_el_vti.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_curv_col_el_aniso.o: forward/sv_curv_col_el_aniso.c
	${CC} -c -o $@ $(CFLAGS) $<
sv_curv_col_ac_iso.o:   forward/sv_curv_col_ac_iso.c
	${CC} -c -o $@ $(CFLAGS) $<

main_curv_col_2d.o: forward/main_curv_col_2d.c
	${CC} -c -o $@ $(CFLAGS) $<

cleanexe:
	rm -f main_curv_col_2d

cleanobj:
	rm -f *.o
cleanall: cleanexe cleanobj
	echo "clean all"
distclean: cleanexe cleanobj
	echo "clean all"
