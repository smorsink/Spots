#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
#base="/Users/kitung/Desktop/thesis_material"
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-25"
#pwd
make spot
times
out_dir="$exe_dir/BBtests"

# integers
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=16     # phase bins
numbands=301   # energy bands
spectraltype=0 # 
#beaming=10     # McPhac
beaming=0      # Blackbody no beaming
spotmodel=0    # circular in the static frame, no gamma
inst_res=0
attenuation=0  # 

# doubles -- using "e" notation
spin=173.6 # in Hz
mass=1.44 # in Msun
radius=13 # in km
inclination=42 # in degrees
emission=56 # in degrees
rho=0.013  # in radians
rho2=0.36  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.231139 # in keV (Surface temperature)
deltatheta=7.78
#deltatheta=0
distance=0.1563 # distance in kpc
elo=0.095
ehi=3.105
obstime=1000000.0      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)
dsback=0.3
agnback=0.4
phase_2=0.5625 # second spot is 0.5625 cycles ahead of first spot
#phase_2=0.565
temp_2=0.0577846
nh=5.0     # in multiples of 4e19 or 1.8e20

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi



# TEST 1: First spot
out_file="$out_dir/one_blackbody.txt"
numtheta=1
spectraltype=0
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot1.txt"
times


# TEST 2: second spot - No Phase Shift 
out_file="$out_dir/two_blackbody.txt"
numtheta=30
spectraltype=0
emission2=131.780293
phase_2=0
rho2=0.36
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission2" -E "$deltatheta" -l "$phase_2" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho2" -T "$temp_2" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot2_30.txt"
times


# TEST 3: First spot - Hopf Beaming
out_file="$out_dir/one_hopf.txt"
numtheta=1
spectraltype=0
beaming=7  #Hopf Beaming
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot1hopf.txt"
times


# TEST 4: second spot - No Phase Shift - Hopf Beaming
out_file="$out_dir/two_hopf.txt"
numtheta=30
spectraltype=0
beaming=7  #Hopf Beaming
emission2=131.780293
phase_2=0
rho2=0.36
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission2" -E "$deltatheta" -l "$phase_2" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho2" -T "$temp_2" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot2_30hopf.txt"
times


# TEST 5: First spot
out_file="$out_dir/one_blackbodytab.txt"
numtheta=1
spectraltype=0
beaming=12
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot1tab.txt"
times

# TEST 6: First spot - Hopf Beaming - Tabulated
out_file="$out_dir/one_hopftab.txt"
numtheta=1
spectraltype=0
beaming=13
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot1hopftab.txt"
times

# TEST 6: First spot - Hopf Beaming - 2D Tabulated
out_file="$out_dir/one_hopftab.txt"
numtheta=1
spectraltype=0
beaming=14
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_dir/spot1hopftab2d.txt"
times