#!/bin/bash
# Author: crwang
# This script can find the KPOINTS number in */OUTCAR
# To use it: ./irreducible_all.sh

for i in *
do
  if [ -e $i/OUTCAR ]
  then
    echo -e $i,$(grep irreducible $i/OUTCAR | awk '{print $2}')
  fi
done
