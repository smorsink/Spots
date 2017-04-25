#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/kitung/Desktop/thesis_material"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/input"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=18     # number of theta bins; for a small spot, only need one
NS_model=3      # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=32      # phase bins
numbands=15     # energy bands
spectraltype=3  # Integrated
elo=0.1
ehi=3.1
beaming=5       # NSXH
spotmodel=0     # circular in the static frame, no gamma
obstime=1000000.0      #length of observation (in seconds)
back=0.001    #constant background (same in all energy bands)

# doubles -- using "e" notation
spin=300 # in Hz
mass=1.48 # in Msun
radius=11.8 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.09669 # in keV (Surface temperature)
distance=0.3 # distance in kpc

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1: Cole's McPHAC, 1.1e6K, 90/90, integrated
out_file="$out_dir/nsxh_apr_19_inte.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times
