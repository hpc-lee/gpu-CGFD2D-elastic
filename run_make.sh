#!/bin/bash

#-- SM code of A100
export SMCODE=sm_80

#-- unknown host
export CUDAHOME=/data3/lihl/software/cuda-11.5
export GNUHOME=/data3/lihl/software/gcc-10.3.0-compile
export NETCDF=/data3/lihl/software/disable-netcdf-4.4.1

echo
echo "start to make ..."
make -j 
