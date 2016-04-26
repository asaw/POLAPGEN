#!/bin/bash
#PBS -l ncpus=1,walltime=100:00:00
#PBS -l mem=1000MB


/home/users/asawigr/reef/R/R-2.13.0/bin/Rscript -e "id_wave_ref <- \"$id_wave\"; id_ch <- \"$id_ch\"; source(\"/home/users/asawigr/reef/zad15/3_wysiew_liscie/1_ref_cow_grid.r\")" >outtest3.txt




