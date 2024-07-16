# Usage: ./check.sh
rm -rf ./tmp

ntot=0
nfinished=0

for i in */
do
    ((ntot++))
    echo '-------------------' >> ./tmp
    if [ -e "${i}OSZICAR" ]; then
        ((nfinished++))
        step=$(grep "=" "${i}OSZICAR" | tail -1)
        conv=$(grep ":" "${i}OSZICAR" | tail -1 | awk '{print $2,$4}')
        gprra=$(grep "required" "${i}OUTCAR")
        echo "$i" "$step" >> ./tmp
        echo "$i" "$conv" >> ./tmp
        echo "$i" "$gprra" >> ./tmp
    else
        echo "$i" '! ! !' >> ./tmp
    fi
done

echo "Total directories: $ntot"
echo "Unfinished directories: $((ntot - nfinished))"

