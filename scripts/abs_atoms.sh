for i in *.vasp
do
  name=$(basename $i .vasp)
  mkdir $name
  python abs_atoms.py $i POSCAR
  mv POSCAR $name/ 
  rm $i
done
