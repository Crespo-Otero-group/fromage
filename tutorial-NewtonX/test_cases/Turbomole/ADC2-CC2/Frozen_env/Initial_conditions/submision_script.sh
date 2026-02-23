#!/bin/sh
#
for i in {1..30}
 do
   cp script.sh I$i/
   cd I$i/
   sed -i "s/NUMB/${i}/g" script.sh
   qsub script.sh
   cd ../
 done
