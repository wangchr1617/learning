# Usage: ./CONTCAR.sh
mkdir ../CONTCAR
for i in */
do
  name=$(basename $i /) 
  cp $name/CONTCAR ../CONTCAR/$name.vasp
done
