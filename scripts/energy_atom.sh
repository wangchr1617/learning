# Usage: ./energy_atom.sh 15
ato=$1
for i in *
do
  if [ -e $i/OUTCAR ]
  then
    ene=$(grep " without" $i/OUTCAR | tail -n 1 | awk '{print $7}')
    ave=$(echo "scale=4; $ene / $ato" | bc)
    echo -e $i,$ave
  fi
done
