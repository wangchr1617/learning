root_path=`pwd`
for cij in `ls -F | grep /$`
do
  cd ${root_path}/$cij
  for s in strain_*
  do
    cd ${root_path}/$cij/$s
    echo `pwd`
  ln -s ${root_path}/opt_ion.py  ./
  ln -s ${root_path}/nep.txt ./
  python opt_ion.py > OUTCAR
  done
done
