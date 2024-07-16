#!/bin/bash
# Auther: Wangchr
# ./kpoints_test.sh 4 5 6 7 8 9

for i in $*
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i
  cp INCAR KPOINTS POTCAR POSCAR $i
  s1="$i  $i  $i"
  sed -i "4s/.*/$s1/g" $i/KPOINTS
done
#qsub batchrunvasp.pbs
