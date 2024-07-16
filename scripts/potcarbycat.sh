# Usage: ./potcarbycat.sh B Mo_sv Ti_sv Al
if [ -e POTCAR ]
then
  rm POTCAR
fi
path=~/POTCAR/PBE
for i in "$@"
do
  cat $path/$i/POTCAR >> POTCAR
done
