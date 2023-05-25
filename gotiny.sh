#!/bin/bash

# TESTS The tiny equatorial spot

# Scripts to run NICER code tests -- Sharon's computer settings
times

base="/Users/sharon/code"
exe_dir="$base/Spots"
pwd
make all
times
out_dir="$exe_dir/Test"

# integers
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=32     # phase bins

spectraltype=0 # 
spotmodel=0    # circular in the static frame, no gamma
inst_res=1     # No instrument=0

# doubles -- using "e" notation
spin=300 # in Hz
mass=1.4 # in Msun
radius=12 # in km

#phaseshift=0.098175 # We need a half-bin phaseshift when comparing with Amsterdam (radians)
distance=0.15 # distance in kpc

obstime=1e6      #length of observation (in seconds)


### Energy Bands
# This set of energy bands works well
elo=0.1
ehi=5.0
#numbands=49
numbands=490

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

### Flipped star. Observer and spot are both in Northern Hemisphere
inclination=43 # in degrees

######
###### TEST 1: First spot - Tiny; on equator - Compute with 0 ISM and 0 instrument
######
NN=1
out_file1="$out_dir/hyd-tiny-NH2-128.txt"
#out_file1="$out_dir/test3.txt"
numtheta="$NN"

# One small spot on equator
emission=90 # Degrees
rho=0.0100 # Radians
beaming=11 #Hydrogen

temp=1.0e6 #Kelvin
#temp=8.617333262e-2 # in keV (Surface temperature)
temp=0.08617315 # keV

nh=200 # units of 1e18 cm^-2
# This is 2e20 cm^-2
#nh=0 # No ISM

## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file1" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -Z "$obstime" -A "$nh" -R "$inst_res"


times

numbins=32
out_file2="$out_dir/hyd-tiny-emit.txt"
## RUNNING THE CODE
./emit -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -q "$NS_model" -l "$phaseshift" -n "$numbins" -b "angles1000.txt" -o "$out_file2" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" 




out_file3="$out_dir/hyd-tiny-test$nh.txt"
## RUNNING THE CODE
./detect  -n "$numbins" -o "$out_file3"  -S "$numbands" -u "$elo" -U "$ehi" -A "$nh" -R "$inst_res" -I "$out_file2" -Z "$obstime"




