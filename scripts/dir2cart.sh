for i in *.vasp
do
  echo -e "411 \n3 \n$i" | vaspkit
  rm $i
done

for i in *REV
do
  name=$(basename $i .vasp_REV)
  mv $i $name.vasp
done
