#!/bin/bash
# Author: crwang
# This script can find the fermi energy in */OUTCAR
# To use it: ./fermi_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep E-fermi $i/OUTCAR | tail -n 1 | awk '{print $3}')
  fi
done
