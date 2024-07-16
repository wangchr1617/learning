# Usage: ./sigma_test.sh 0.05 0.1 0.15 0.2

for i in $*
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i

  cp INCAR KPOINTS POTCAR POSCAR $i
  s1="SIGMA = "
  s2="SIGMA = $i #"
  sed -i "s/$s1/$s2/g" $i/INCAR
done
#qsub batchrunvasp.pbs

