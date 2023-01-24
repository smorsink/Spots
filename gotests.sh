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
numbins=32     # phase bins

spectraltype=0 # 
spotmodel=0    # circular in the static frame, no gamma
inst_res=0     # No instrument
attenuation=0  # No ISM

# doubles -- using "e" notation
spin=300 # in Hz
mass=1.4 # in Msun
radius=12 # in km

phaseshift=0.0 # this is used when comparing with data
distance=0.15 # distance in kpc

obstime=1      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)
dsback=0.0
agnback=0.0

### Energy Bands
#elo=0.095
#ehi=3.105

elo=0.1
ehi=3.1
numbands=300


elo=0.1
ehi=4.0
numbands=39

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

### Flipped star. Observer and spot are both in Northern Hemisphere
inclination=43 # in degrees

NN=1

# TEST 1: First spot - Tiny; on equator
spin=300
out_file1="$out_dir/hyd-tiny.txt"
numtheta="$NN"

# One small spot near South Pole
emission=90 # Degrees
rho=0.0100 # Radians
beaming=11

temp=1.0e6 #Kelvin
#temp=8.617333262e-2 # in keV (Surface temperature)
temp=0.08617315 # keV
nh=5

## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model"  -E "$deltatheta" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file1" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" -A "$nh" -R 1






times

