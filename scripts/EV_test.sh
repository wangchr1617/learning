# Usage: ./EV_test.sh

for i in 0.90 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.10
do
  if [ -d $i ]
  then
    rm -r $i
  fi
  mkdir $i
  cp INCAR KPOINTS POSCAR POTCAR *.pbs $i
  sed -i "2s/1.0/$i/g" $i/POSCAR
done

# natoms=`sed -n '7p' POSCAR | awk '{for(i = 1; i <= NF; i++) {sum += $i} {printf("%d", sum)}}'`
# for i in *
# do
#   if [ -e $i/OUTCAR ]
#   then
#     ene=$(grep " without" $i/OUTCAR | tail -n 1 | awk '{print $7}')
#     ave=$(echo "scale=4; $ene / $natoms" | bc)
#     vol=$(grep volume $i/OUTCAR | tail -n 1 | awk '{print $5}')
#     echo -e $vol,$ave
#   fi
# done
