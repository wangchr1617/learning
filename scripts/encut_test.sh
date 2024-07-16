#!/bin/bash
# Auther: Wangchr
# ./encut_test.sh 150 200 250 300

for i in $*
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i

  cp INCAR KPOINTS POTCAR POSCAR $i
  s1="ENCUT = "
  s2="ENCUT = $i #"
  sed -i "s/$s1/$s2/g" $i/INCAR 
done
#qsub batchrunvasp.pbs

