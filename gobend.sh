#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-9"
#pwd
make bend
times
out="angles100.txt"

./bend -l 0.1 -h 0.284 -n 100 -o "$out" 

#./bend -c 0 -r 1 -i "angles50.txt"

times