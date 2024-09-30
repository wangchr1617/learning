# Usage: ./vasp2xyz_all.sh

num=0
for i in *
do
  if [ -e $i/OSZICAR ]
  then
    cd $i
    if [ -e NEP-dataset.xyz ]
    then
        rm NEP-dataset.xyz
    fi
    python vasp2xyz.py
    echo $i
    num=$[$num+1]
    cat NEP-dataset.xyz >> ../NEP-dataset.xyz
    cd $OLDPWD
  fi
done 
echo $num
