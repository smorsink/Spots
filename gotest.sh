#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-11"
#pwd
make spot
times
out_dir="$exe_dir/OBLATE-gamma"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=1     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numbands=1     # energy bands
spectraltype=0 # Blackbody observed at one Energy
elow=0.3
ehi=0.8
beaming=0      # isotropic emission, no beaming
spotmodel=0


# doubles -- using "e" notation
spin=600 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1.0  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.1 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST OS1a = 600 Hz, tiny spot; oblate
numtheta=5 # number of theta bins; for a small spot, only need one
spotmodel=0
NS_model=3
beaming=0
spectraltype=2
out_file="$out_dir/test-5.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt"
times


