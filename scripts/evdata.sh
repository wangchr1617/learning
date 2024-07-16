# Author: crwang
# This script can find the energy/atom in */OUTCAR
# Usage: ./evdata.sh 15

natoms=$1
for i in *
do
  if [ -e $i/OUTCAR ]
  then
    ene=$(grep " without" $i/OUTCAR | tail -n 1 | awk '{print $7}')
    ave=$(echo "scale=4; $ene / $natoms" | bc)
    vol=$(grep volume $i/OUTCAR | tail -n 1 | awk '{print $5}')
    echo -e $vol,$ave
  fi
done
