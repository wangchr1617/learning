# Usage: ./cp2k2vasp.sh
if [ -e cp2k-1.restart ]
then
echo -e "cp2k-1.restart \n100 \n2 \n27 \nCONTCAR \n0 \nq" | Multiwfn_noGUI
fi

