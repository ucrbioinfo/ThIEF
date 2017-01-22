#!/bin/sh
glpsol --lp $1 -o ${1%.lp}.sol
grep 'x' ${1%.lp}.sol  | awk '$4==1 {print substr($2,2,10)}' > ${1%.lp}.csv