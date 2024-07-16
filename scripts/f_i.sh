#!/bin/bash
# Author: crwang
# This script can find the imaginary frequency in */OUTCAR
# To use it: ./f_i.sh
grep 'f/i' */OUTCAR | awk '{print $1 "\t" $2 "\t" $8 "\t" $9 "\t" $10 "\t" $11}'
