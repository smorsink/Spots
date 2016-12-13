#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-12"
#pwd
make spot
times
out_dir="$exe_dir/ML2015_data"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=10     # number of theta bins; for a small spot, only need one
NS_model=2     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=16    # phase bins
numbands=30     # energy bands
spectraltype=2 # Blackbody energy bands
elo=3.5
ehi=12.5
beaming=7      # Hopf function limb-darkening
#spotmodel=0    # circular in the static frame, with gamma
#spotmodel=1    # circular in the rotating frame, with gamma
spotmodel=2    # circular in the static frame, no gamma
obstime=1040500.0      #length of observation (in seconds)
back=0.0705    #constant background (same in all energy bands)

# doubles -- using "e" notation
spin=600 # in Hz
mass=1.6 # in Msun
radius=11.8 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.4363  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=2.0 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
#distance=6.1713552e18 # 200 pc in meters
distance=10.0 # distance in kpc

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST OS1a = 600 Hz, tiny spot; oblate

out_file="$out_dir/testa.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
#times

# TEST OS1b = 600 Hz, big spot
NS_model=1
spotmodel=0
numbins=128
numtheta=100
rho=1.0
out_file="$out_dir/os1b.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

