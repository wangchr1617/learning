#!/bin/bash
# Auther: Wangchr
# ./encut_test.sh 150 200 250 300

natom=`sed -n '7p' POSCAR | awk '{for(i = 1; i <= NF; i++) {sum += $i} {printf("%d", sum)}}'`
echo $natom

for i in $*
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i

  cp INCAR POTCAR POSCAR $i
  s1="ENCUT = "
  s2="ENCUT = $i #"
  sed -i "s/$s1/$s2/g" $i/INCAR 
done
