#!/bin/bash

#-- SM code of A100
export SMCODE=sm_80

export GNUHOME=/usr
export CUDAHOME=/usr/local/cuda-11.8
export NETCDF=/data/apps/NetCDF/disable-netcdf-4.8.1

echo
echo "start to make ..."
make -j 
