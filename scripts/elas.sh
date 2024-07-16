i='./'
for i in */
do
  cd $i
  if [ -e ./tmp ]
  then
    rm ./tmp
  fi
  echo -e "200" | vaspkit >> ./tmp
  loc=`pwd`
  echo $(basename $loc /) >> ../elas_cons
  grep -A3 "Stiffness Tensor" ./tmp >> ../elas_cons 
  echo "--------------------------" >> ../elas_cons
  cd ..
done
