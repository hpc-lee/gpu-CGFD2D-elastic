#******************************************************************************#
#*                                                                            *#
#*                               CGFD2D-Wave                                  *#
#*                           =========================                        *#
#*                                                                            *#
#*   Curvilinear Grid Finite-Difference 2D Elastic/Acoustic Wave modeling     *#
#*                                                                            *#
#*     A signed non-commercial agreement is required to use this program.     *#
#*             Free for non-commercial academic research ONLY.                *#
#*       This program is distributed WITHOUT ANY WARRANTY whatsoever.         *#
#*       Do not redistribute this program without written permission.         *#
#*                                                                            *#
#* Written by:                                                                *#
#*   Wei ZHANG      (zhangwei@sustech.edu.cn)                                 *#
#*                                                                            *#
#* Revision History:                                                          *#
#*   2021/08/24     original version modified from CGFD3D-Wave                *#
#******************************************************************************#

1 Introduction

  This program simulates elastic/asoutic iso/aniso wave propagation in 2D complex
  media with surface topography by using curvilinear-grid finite-difference method,
  and with flat surface by using staggered-grid finite-difference sheme.
  Details about the algorithm can be found in Zhang and Chen (2006).

2 Direcotry layout

    Makefile 
       makefile for this package

    example/
       run script to demonstrate the usage
       matlab scriptes to show the result
       #need snctools: http://mexcdf.sourceforge.net/

    forward/
       source codes

    lib/
       lib source codes

    media/
       codes for discreting media, written by Luqian Jiang 

    README
       this file

3 Installation

    3.1 You need to install a c compiler, NetCDF

    3.2 Compiling:
        
        edit Makefile
            there are three executables:
                main_curv_col_2d : curvilinear grid, macdrp schemem, for el/ac/aniso media
                main_cart_col_el_2d: cartesian grid, macdrp schemem, only el iso media
                main_cart_stg_2d: cartesian grid, staggered scheme, for el/ac + iso media

        make all

4 Usage
  
  1) cd example/
  2) edit run_test.sh
  3) run run_test.sh
  4) use .m to show results

# vim:ft=conf:ts=4:sw=4:nu:et:ai:
