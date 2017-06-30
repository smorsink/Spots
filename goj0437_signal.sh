#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/kitung/Desktop/thesis_material"
#base="/Users/sharon/code/Albert"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/team_hw"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=10    # number of theta bins; for a small spot, only need one
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=16     # phase bins
numbands=301   # energy bands
spectraltype=0 # 
beaming=11     # NSXHnew
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
temp_2=0.0577846
nh=5.0     # in multiples of 4e19 or 1.8e20

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1: NSXH full table
out_file="$out_dir/jun30_nsxhnew_obl_j0437_mono.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times

# TEST 1b: NSXH full table
out_file="$out_dir/jun30_nsxhnew_obl_j0437_inte.txt"
spectraltype=3
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times

# TEST 2: Cole's McPHAC
out_file="$out_dir/jun30_mcphacc_obl_j0437_mono.txt"
spectraltype=0
beaming=10     # McPhac
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -2 -B "$phase_2" -C "$temp_2" -d "$rho2" -O "$out_dir/bothspots30.txt"
times

# TEST 2b: Cole's McPHAC
out_file="$out_dir/jun30_mcphacc_obl_j0437_inte.txt"
spectraltype=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times

# TEST 3: NSX Helium full table
out_file="$out_dir/jun30_nsxhenew_obl_j0437_mono.txt"
spectraltype=0
beaming=15     # NSX Helium new full table
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times

# TEST 3b: NSX Helium full table
out_file="$out_dir/jun30_nsxhenew_obl_j0437_inte.txt"
spectraltype=3
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times
