#!/bin/sh
#
for i in {1..10}
 do
   cp script.sh TRAJ$i/
   cp fromage.in TRAJ$i/JOB_AD/
   cd TRAJ$i/
   sed -i "s/NUMB/${i}/g" script.sh
   qsub script.sh
   cd ../
 done
