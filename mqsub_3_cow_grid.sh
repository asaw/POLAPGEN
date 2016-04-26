#!/bin/bash
for (( i = 1; $i <= 105; i++ ))
do
  for j in 280 330
  do
     id_variety=`printf %03d $i`
      qsub -v id_variety=$i,id_wave=$j /home/users/asawigr/reef/zad15/3_wysiew_liscie/run_3_cow_grid.sh
  done
done
