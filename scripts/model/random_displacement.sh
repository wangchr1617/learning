for i in {1..500}
do
  random_int=$(awk -v min=1 -v max=400 'BEGIN{srand(); print int(min+rand()*(max-min+1))}')
  random_decimal=$(awk "BEGIN { printf \"%.4f\", $random_int/10000 }")
  a=$(awk "BEGIN { printf \"%.4f\", $random_decimal + 0.01 }")
  echo -e "407 \n1 \n1 \nGe \n$a" | vaspkit
  mv POSCAR_REV ${a}_${i}.vasp
done
