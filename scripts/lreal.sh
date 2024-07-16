for a in */
do
natom=`sed -n '7p' $a/POSCAR | awk '{for(i = 1; i <= NF; i++) {sum += $i} {printf("%d", sum)}}'`
echo $a, $natom
if test $(($natom)) -gt 100; then
sed -i 's/LREAL = F/LREAL = A/' $a/INCAR
else
sed -i 's/LREAL = A/LREAL = F/' $a/INCAR
fi
done

