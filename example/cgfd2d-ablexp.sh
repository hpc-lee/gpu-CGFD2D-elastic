#!/bin/bash

#set -x
set -e

date

#-- program related dir
EXEC_WAVE=`pwd`/../main_curv_col_2d
echo "EXEC_WAVE=$EXEC_WAVE"

#-- input dir
INPUTDIR=`pwd`

#-- output and conf
PROJDIR=`pwd`/../project
PAR_FILE=${PROJDIR}/test.json
GRID_DIR=${PROJDIR}/output
MEDIA_DIR=${PROJDIR}/output
SOURCE_DIR=${PROJDIR}/output
OUTPUT_DIR=${PROJDIR}/output

rm -rf ${PROJDIR}

#-- create dir
mkdir -p $PROJDIR
mkdir -p $OUTPUT_DIR
mkdir -p $GRID_DIR
mkdir -p $MEDIA_DIR

#----------------------------------------------------------------------
#-- grid and mpi configurations
#----------------------------------------------------------------------

#-- total x grid points
NX=801
#-- total z grid points
NZ=401
#----------------------------------------------------------------------
#-- create main conf
#----------------------------------------------------------------------
cat << ieof > $PAR_FILE
{
  "number_of_total_grid_points_x" : ${NX},
  "number_of_total_grid_points_z" : ${NZ},

  "#size_of_time_step" : 0.00005,
  "#number_of_time_steps" : 10000,
  "time_window_length" : 0.5,
  "check_stability" : 1,

  "boundary_x_left" : {
      "ablexp" : {
          "number_of_layers" : 50,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_x_right" : {
      "ablexp" : {
          "number_of_layers" : 50,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_z_bottom" : {
      "ablexp" : {
          "number_of_layers" : 50,
          "ref_vel"  : 7000.0
          }
      },
  "boundary_z_top" : {
      "free" : "timg"
      },

  "grid_generation_method" : {
      "import" : "$INPUTDIR/grid2",
      "#cartesian" : {
        "origin"  : [0.0, -400.0 ],
        "inteval" : [ 1.0, 1.0 ]
      }
  },
  "is_export_grid" : 1,
  "grid_export_dir"   : "$GRID_DIR",

  "metric_calculation_method" : {
      "#import" : "$GRID_DIR",
      "calculate" : 1
  },
  "is_export_metric" : 1,

  "medium" : {
      "type" : "elastic_iso",
      "#input_way" : "infile_layer",
      "#input_way" : "binfile",
      "input_way" : "code",
      "#binfile" : {
        "size"    : [1101, 1447, 1252],
        "spacing" : [-10, 10, 10],
        "origin"  : [0.0,0.0,0.0],
        "dim1" : "z",
        "dim2" : "x",
        "dim3" : "y",
        "Vp" : "$INPUTDIR/prep_medium/seam_Vp.bin",
        "Vs" : "$INPUTDIR/prep_medium/seam_Vs.bin",
        "rho" : "$INPUTDIR/prep_medium/seam_rho.bin"
      },
      "code" : "func_name_here",
      "#import" : "$MEDIA_DIR",
      "#infile_layer" : "$INPUTDIR/prep_medium/basin_el_iso.md3lay",
      "#infile_grid" : "$INPUTDIR/prep_medium/topolay_el_iso.md3grd",
      "#equivalent_medium_method" : "loc",
      "#equivalent_medium_method" : "har"
  },

  "is_export_media" : 1,
  "media_export_dir"  : "$MEDIA_DIR",

  "#visco_config" : {
      "type" : "graves_Qs",
      "Qs_freq" : 1.0
  },

  "in_source_file" : "$INPUTDIR/test_source.src",
  "is_export_source" : 1,
  "source_export_dir"  : "$SOURCE_DIR",

  "output_dir" : "$OUTPUT_DIR",

  "in_station_file" : "$INPUTDIR/station.list",

  "#receiver_line" : [
    {
      "name" : "line_x_1",
      "grid_index_start"    : [  50, 59 ],
      "grid_index_incre"    : [  5,  0 ],
      "grid_index_count"    : 10
    },
    {
      "name" : "line_z_1",
      "grid_index_start"    : [ 200, 39 ],
      "grid_index_incre"    : [  0,  2 ],
      "grid_index_count"    : 10
    } 
  ],

  "#snapshot" : [
    {
      "name" : "volume_vel",
      "grid_index_start" : [ 0,   0 ],
      "grid_index_count" : [ ${NX}, ${NZ} ],
      "grid_index_incre" : [  1,  1 ],
      "time_index_start" : 0,
      "time_index_incre" : 1,
      "save_velocity" : 1,
      "save_stress"   : 0,
      "save_strain"   : 0
    }
  ],

  "check_nan_every_nummber_of_steps" : 0,
  "output_all" : 0 
}
ieof

echo "+ created $PAR_FILE"

#-------------------------------------------------------------------------------
#-- Performce simulation
#-------------------------------------------------------------------------------
#
#-- gen run script
cat << ieof > ${PROJDIR}/cgfd_sim.sh
#!/bin/bash

set -e

printf "\nStart simualtion ...\n";
time $EXEC_WAVE $PAR_FILE 100 3 2>&1 |tee log
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

ieof

#-------------------------------------------------------------------------------
#-- start run
#-------------------------------------------------------------------------------

chmod 755 ${PROJDIR}/cgfd_sim.sh
${PROJDIR}/cgfd_sim.sh
if [ $? -ne 0 ]; then
    printf "\nSimulation fail! stop!\n"
    exit 1
fi

date

# vim:ts=4:sw=4:nu:et:ai:
