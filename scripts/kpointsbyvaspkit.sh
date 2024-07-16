# Usage: ./kpointsbyvaspkit.sh
for i in */
do
  if [ -e $i/INCAR ]
  then
    cd $i
    echo -e "102 \n2 \n0.04" | vaspkit
    cd $OLDPWD
  fi
done
