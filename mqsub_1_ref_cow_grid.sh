#!/bin/bash
  for j in 280 330
  do
    for (( k = 1; $k <= 105; k++ ))
    do
     id_ch=`printf %03d $k`
     qsub -v id_wave=$j,id_ch=$k /home/users/asawigr/reef/zad15/3_wysiew_liscie/run_1_ref_cow_grid.sh
    done
  done