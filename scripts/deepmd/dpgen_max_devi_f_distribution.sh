if [ -e fdr.dat ]
then
  rm fdr.dat
fi
find ./ -name "model_devi.out" >> fdr.dat
for t in 300 450 600 750
do
  if [ -e ${t}K ]
  then
    rm ${t}K
  fi
  for i in `cat fdr.dat`
  do
    cat $i | awk '{print $5}' >> ${t}K
  done
  echo "Finished collection of temperature $t."
  sed -ie '/avg/d' ${t}K
done
paste *K > Max_Devi_F 
