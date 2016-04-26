#!/bin/bash
#PBS -l ncpus=1,walltime=24:00:00
#PBS -l mem=500MB

/home/users/asawigr/reef/R/R-2.13.0/bin/Rscript -e "id_wave <- \"$id_wave\"; source(\"/home/users/asawigr/reef/zad15/3_wysiew_liscie/2_ref_cow_opt.r\")" >outtest4.txt




