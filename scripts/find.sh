# Usage: ./find.sh CHGCAR WAVECAR
for i in "$@"
do
  find ./ -name "$i" >> name.sh
done

sed -i 's/^/rm /' name.sh
# sed -i 's/^/cat /' name.sh
# sed -i 's/$/ >> ~\/unary.xyz/' name.sh

chmod u+x name.sh
./name.sh
rm name.sh 
