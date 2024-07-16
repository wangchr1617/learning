mkdir $1
for i in *.vasp
do
  cp $i $1/$1-$i
  sed -i "2s/1.0/$1/g" $1/$1-$i
done
