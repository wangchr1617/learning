# Usage: 
for i in */
do
  if [ -e $i/INCAR ]
  then
    cd $i
    echo -e "102 \n2 \n0.04" | vaspkit # KPOINTS
    cd $OLDPWD
  fi
done
