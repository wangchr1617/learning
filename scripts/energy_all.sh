#!/bin/bash
# Author: crwang
# This script can find the energy in */OUTCAR
# To use it: ./energy_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep " without" $i/OUTCAR | tail -n 1 | awk '{print $7}')
  fi
done
