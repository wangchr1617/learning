#!/bin/bash
# Author: crwang
# This script can submit tasks of all folders
# To use it: ./qsub_all.sh

for i in *
do
  if [ -e $i/INCAR ]
  then
    cd $i
    qsub runvasp.pbs
    cd $OLDPWD
  fi
done
