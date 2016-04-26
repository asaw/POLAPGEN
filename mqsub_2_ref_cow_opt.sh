#!/bin/bash
  for j in 280 330
  do
     qsub -v id_wave=$j /home/users/asawigr/reef/zad15/3_wysiew_liscie/run_2_ref_cow_opt.sh
  done