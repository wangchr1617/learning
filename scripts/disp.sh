# Usage: ./disp.sh
for i in POSCAR-0*
do
  name=$(echo "$i" | grep -oE '[0-9]+') 
  echo $name
  mkdir $name
  mv $i $name/POSCAR
done
