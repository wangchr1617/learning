rm -rf ./tmp
for i in */
do
echo '-------------------' >> ./tmp
if [ -e $i/output ]; then
gptim=`grep "Time used =" $i/output`
gpfin=`grep "Finished running" $i/output`
gpdum=`grep = $i/dump.xyz | wc -l `
echo $i $gptim >> ./tmp
echo $i $gpfin >> ./tmp
echo $i $gpdum >> ./tmp
else
echo $i '! ! !' >> ./tmp
fi
done

