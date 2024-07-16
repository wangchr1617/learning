#!/bin/bash
# Author: crwang
# This script can find the frequency in */OUTCAR and sum them. The output is in eV.
# To use it: ./fre_sum.sh

hv_sum=$(grep "f  =" OUTCAR | awk '{print $10}' | paste -sd+ | bc)
echo "scale = 6; $hv_sum/2000" | bc
