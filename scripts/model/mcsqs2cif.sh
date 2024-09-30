# Usage: ./mcsqs2cif.sh
for i in bestsqs*.out
do
  name=$(basename $i .out)
  str2cif < $i > $name.cif
done

