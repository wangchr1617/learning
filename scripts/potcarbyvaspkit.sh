for i in *
do
  if [ -e $i/INCAR ]
  then
    cd $i
    echo -e "103 \n" | vaspkit
    cd $OLDPWD
  fi
done
