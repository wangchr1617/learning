#!/bin/bash
# Author: crwang
# This script can find the freq in OUTCAR
# To use it: ./freq.sh

echo -e $(grep "f  =" OUTCAR | awk '{print $10}')
