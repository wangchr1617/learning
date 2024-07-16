#!/bin/bash
# Author: crwang
# This script can find the eentro in */OUTCAR
# To use it: ./eentro_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep "EENTRO" $i/OUTCAR | tail -n 1 | awk '{print $5}')
  fi
done
