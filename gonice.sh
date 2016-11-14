#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-6"
#pwd
make spot
times
out_dir="$exe_dir/NICER"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numbands=1     # energy bands
spectraltype=0 # Blackbody observed at one Energy
beaming=0      # Blackbody, no beaming


# doubles -- using "e" notation
spin=1 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1e-2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.35 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1 = 1 Hz, tiny spot
numtheta=1 # number of theta bins; for a small spot, only need one
out_file="$out_dir/test1.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands"
#times

# TEST 2 = 1 Hz, big spot
numtheta=50
rho=1
out_file="$out_dir/test2.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands"
#times

# TEST 3 = 200 Hz, tiny spot
numtheta=1 # number of theta bins; for a small spot, only need one
#numphi=1
numbins=128
spin=200
rho=1e-2
out_file="$out_dir/test3.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands"
#times

# TEST 4 = 200 Hz, big spot
numtheta=50 # number of theta bins; for a small spot, only need one
numphi=16
spin=200
rho=1
out_file="$out_dir/test4.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands"
times

# TEST 5 = 400 Hz, big spot
numbins=128
numtheta=50 # number of theta bins; for a small spot, only need one
spin=400
rho=1
inclination=30 # in degrees
emission=60 # in degrees
beaming=0
out_file="$out_dir/test5.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands"
times

# TEST 6 = 400 Hz, big spot

numbins=128
numtheta=50 # number of theta bins; for a small spot, only need one
numphi=20
spin=400
rho=1
inclination=80 # in degrees
emission=20 # in degrees
#beaming=4   # Helium
out_file="$out_dir/test6.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -g "$beaming"

times
