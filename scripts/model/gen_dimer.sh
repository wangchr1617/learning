dimer="echo GeGe"
nebmake.pl ini mid 9 

for i in {00..10}
do
  name=$(expr "scale=1;1.0 + $i/10"|bc)
  fname=`$dimer`_$name
  mv $i $fname
  cp INCAR POTCAR $fname
done

nebmake.pl mid fin 24

rm -r 00/
for j in {01..25}
do 
  name=$(expr "scale=1;2.0 + $j/5"|bc)
  fname=`$dimer`_$name
  mv $j $fname
  cp INCAR POTCAR $fname
done

