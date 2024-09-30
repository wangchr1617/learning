num=0
for i in *
do
  if [ -e $i/OUTCAR ]
  then
    cd $i
    ./singleFrame-outcars2nep-exyz.sh ./
    num=$[$num+1]
    cat NEPdataset/NEP-dataset.xyz >> ../NEP-dataset.xyz
    cd $OLDPWD
  fi
done 
echo $num
