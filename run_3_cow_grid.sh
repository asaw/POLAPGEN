#!/bin/bash
#PBS -l ncpus=1,walltime=1000:00:00
#PBS -l mem=5000MB

/home/users/asawigr/reef/R/R-2.13.0/bin/Rscript -e "id_wave <- \"$id_wave\"; id_variety <- \"`printf %03d $id_variety`\"; source(\"/home/users/asawigr/reef/zad15/3_wysiew_liscie/3_cow_grid.r\")" >outtest5.txt




