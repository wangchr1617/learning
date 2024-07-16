#!/bin/bash
# Author: crwang
# This script can find the volume in */OUTCAR
# To use it: ./volume_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep volume $i/OUTCAR | tail -n 1 | awk '{print $5}')
  fi
done
