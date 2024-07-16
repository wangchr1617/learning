#!/bin/bash
# Author: crwang
# This script can find the energy in */OUTCAR
# To use it: ./energy_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep "TOTEN" $i/OUTCAR | tail -n 1 | awk '{print $5}')
  fi
done
