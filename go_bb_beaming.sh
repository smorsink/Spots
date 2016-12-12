#!/bin/bash

times
base="/home/kitung"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/bb_3_beaming"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=5     # number of theta bins; for a small spot, only need one
NS_model=2     # 1 (oblate 2014 A&M), 2 (oblate 2007 MLCB), or 3 (spherical);
numbins=16    # phase bins
numbands=5     # energy bands
spectraltype=2 # Integrated Blackbody
elow=0.3	   # Lowest energy boundary
ehi=1.8		   # Highest energy boundary
beaming=0      # isotropic emission, no beaming
spotmodel=0	   # circular spot


# doubles -- using "e" notation
spin=600 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1.0  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=2 # in keV (Surface temperature)
distance=3.08567758e20 # 10 kpc in meters
#distance=6.1713552e18 # 200 pc in meters


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1
numtheta=5 # number of theta bins; for a small spot, only need one
spotmodel=0
NS_model=3
beaming=0
spectraltype=2
out_file="$out_dir/bb_only.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt"
times

# Test 1a
out_file="$out_dir/bb_only_2.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt" -2
times

# TEST 2
beaming=1
out_file="$out_dir/bb_chan.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt"
times

# TEST 2a
out_file="$out_dir/bb_chan_2.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt" -2
times

# TEST 3
beaming=7
out_file="$out_dir/bb_hopf.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt"
times

# TEST 3a
out_file="$out_dir/bb_hopf_2.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt" -2
times

# TEST 4
out_file="$out_dir/hyd.txt"
spectraltype=3
beaming=3
temp=0.1
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -P "$spotmodel" -b "angles100.txt" -2
times

