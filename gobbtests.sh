#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times

base="/Users/sharon/code"
exe_dir="$base/Spots"
pwd
make spot
times
out_dir="$exe_dir/Test"

# integers
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=64     # phase bins
numbands=301   # energy bands
spectraltype=0 # 
#beaming=10     # McPhac
beaming=0      # Blackbody no beaming
spotmodel=0    # circular in the static frame, no gamma
inst_res=0     # No instrument
attenuation=0  # No ISM

# doubles -- using "e" notation
spin=300 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=43 # in degrees
emission=90 # in degrees
rho=0.01  # in radians
rho2=0.36  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=8.617333262e-2 # in keV (Surface temperature)
deltatheta=7.78
#deltatheta=0
distance=0.15 # distance in kpc
elo=0.095
ehi=3.105
obstime=1000000.0      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)
dsback=0.0
agnback=0.0
phase_2=0.5625 # second spot is 0.5625 cycles ahead of first spot
#phase_2=0.565
temp_2=0.0577846
nh=5.0     # in multiples of 4e19 or 1.8e20

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

### Flipped star. Observer and spot are both in Northern Hemisphere
distance=3.2407557e-21 # Remove Distance
distance=1.0
obstime=1.0



NN=128

# TEST 1: First spot - Hydrogen; Medium, on equator
out_file1="$out_dir/hyd-meda-$NN.txt"
numtheta="$NN"
spectraltype=0
# One small spot near South Pole
emission=90 # Degrees
rho=0.5 # Radians
beaming=11
phaseshift=0.2
temp=1e6


## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file1" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 





# TEST 2: Second spot - Hydrogen; Medium, on equator
out_file2a="$out_dir/hyd-med2a-$NN.txt"
numtheta="$NN"
spectraltype=0
# One small spot near South Pole
emission=90 # Degrees
rho=0.25 # Radians
beaming=11
phaseshift=0.0 # Radians


## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file2a" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 

# TEST 3: Second spot - Hydrogen; Medium, on equator
out_file2b="$out_dir/hyd-med2b-$NN.txt"
numtheta="$NN"
spectraltype=0
# One small spot near South Pole
emission=90 # Degrees
rho=0.25 # Radians
beaming=11
phaseshift=0.2


## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file2b" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 



out_file3a="$out_dir/hyd-a-sum$NN.txt"
./addcurves -t 64 -e 300 -i "$out_file1" -j "$out_file2a" -o "$out_file3a"

out_file3b="$out_dir/hyd-b-sum$NN.txt"
./addcurves -t 64 -e 300 -i "$out_file1" -j "$out_file2b" -o "$out_file3b"


times

