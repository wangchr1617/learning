for i in "$@"
do
  for j in */
  do
    ln -s ../$i $j/$i
  done
done

