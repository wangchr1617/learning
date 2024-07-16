#!/bin/bash
# Author: crwang
# This script can print calculation time of different * from */OUTCAR
# To use it: ./time_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo $i,$(grep User $i/OUTCAR | awk '{print $4}')
  fi
done
