# Usage: ./POSCAR.sh GeTe/
cd $1
for i in *
do
  name=$(basename $i .vasp) 
  mkdir $name
  mv $i $name/POSCAR
done
cd ..
