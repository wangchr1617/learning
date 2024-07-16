# Author: crwang
# This script can change NCORE of different dir
# Usage: ./ncore_test.sh 4 8 12

for i in $*
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i

  cp INCAR KPOINTS POTCAR POSCAR $i
  s1="NCORE = "
  s2="NCORE = $i #"
  sed -i "s/$s1/$s2/g" $i/INCAR
done
#qsub batchrunvasp.pbs

