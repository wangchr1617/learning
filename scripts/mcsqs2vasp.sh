# Usage: ./mcsqs2vasp.sh
for i in bestsqs*.out
do
  name=$(basename $i .out)
  sqs2poscar $i
  mv $i-POSCAR $name.vasp
  sed -i "s/xxx/1.0/g"  $name.vasp
done

