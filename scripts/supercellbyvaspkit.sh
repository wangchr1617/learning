for i in */
do
  if [ -e $i/POSCAR ]
  then
    cd $i
    (echo 401; echo 1; echo 3 3 3) | vaspkit
    mv SC333.vasp CONTCAR
    cd $OLDPWD
  fi
done
