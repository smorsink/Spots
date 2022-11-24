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
distance=3.240779e-20 # Remove Distance


# TEST 1: First spot - Blackbody; Medium, on equator
out_file="$out_dir/bb-med-n8.txt"
numtheta=8
spectraltype=0
beaming=0
# One small spot on South Pole
emission=90 # Degrees
rho=0.5 # Radians
#distance=1.0

## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 


NN=128

# TEST 2: First spot - Hydrogen; Medium, on equator
out_file="$out_dir/hyd-med-$NN.txt"
numtheta="$NN"
spectraltype=0
# One small spot near South Pole
emission=90 # Degrees
rho=0.5 # Radians
beaming=11


## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 


# TEST 3: Second spot - Hydrogen; Medium, on equator
out_file2="$out_dir/hyd-med2-$NN.txt"
numtheta="$NN"
spectraltype=0
# One small spot near South Pole
emission=90 # Degrees
rho=0.25 # Radians
beaming=11
phaseshift=0.0318309


## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file2" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" 


out_file3="$out_dir/hyd-sum$NN.txt"

./addcurves -t 64 -e 300 -i "$out_file" -j "$out_file2" -o "$out_file3"


times

