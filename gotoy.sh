#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/kitung/Desktop/thesis_material"
#base="/Users/sharon/code/Albert"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/toy_model"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=10    # number of theta bins; for a small spot, only need one
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=32     # phase bins
numbands=301   # energy bands
spectraltype=3 # integrated
beaming=11     # NSXHnew
spotmodel=0    # circular in the static frame, no gamma
inst_res=1
attenuation=0  # 

# doubles -- using "e" notation
spin=200 # in Hz
mass=1.4 # in Msun
radius=10 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.1  # in radians
#rho2=0.36  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.0967 # in keV (Surface temperature)
#deltatheta=7.78
#deltatheta=0
distance=0.1563 # distance in kpc
elo=0.095
ehi=3.105
obstime=1000000.0      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)
#dsback=0.3
#agnback=0.4
#phase_2=0.5625 # second spot is 0.5625 cycles ahead of first spot
#temp_2=0.0577846
nh=5.0     # in multiples of 4e19

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 0: Default
out_file="$out_dir/jul6_nsxhnew_default.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 1: High Spin
out_file="$out_dir/jul6_nsxhnew_high_spin.txt"
spin=600
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 2: Low Temperature
out_file="$out_dir/jul6_nsxhnew_low_temp.txt"
spin=200
temp=0.0544 #10^5.8 K
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 2a: High Temperature
out_file="$out_dir/jul6_nsxhnew_high_temp.txt"
temp=0.1719 #10^6.3 K
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 3: Higher compactness/M
out_file="$out_dir/jul6_nsxhnew_high_M.txt"
temp=0.0967
mass=1.8 # in Msun
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 4: Larger Radius
out_file="$out_dir/jul6_nsxhnew_high_R.txt"
mass=1.4 # in Msun
rho=0.08
radius=12.5 # in km
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 4: Smaller Radius
out_file="$out_dir/jul6_nsxhnew_low_R.txt"
radius=8 # in km
rho=0.125
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 5: Second Spot with lower temperature
out_file="$out_dir/jul6_nsxhnew_2_spot.txt"
radius=10 # in km
rho=0.1
rho2=0.2  # in radians
deltatheta=0
phase_2=0.5 # second spot is 0.5 cycles ahead of first spot
temp_2=0.0544 #10^5.8 K
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2 -B "$phase_2" -C "$temp_2" -d "$rho2"
times

# TEST 6: Alternative Spot Position
out_file="$out_dir/jul6_nsxhnew_9070.txt"
inclination=90 # in degrees
emission=70 # in degrees
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times 

# TEST 6a: Alternative Spot Position
out_file="$out_dir/jul6_nsxhnew_3090.txt"
inclination=30 # in degrees
emission=90 # in degrees
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times 

# TEST 6b: Alternative Spot Position
out_file="$out_dir/jul6_nsxhnew_3070.txt"
inclination=30 # in degrees
emission=70 # in degrees
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times 

# TEST 7: With normal ISM
out_file="$out_dir/jul6_nsxhnew_ism_2e20.txt"
inclination=90 # in degrees
emission=90 # in degrees
attenuation=5
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuation" -A "$nh"
times 

# TEST 7a: With less ISM
out_file="$out_dir/jul6_nsxhnew_ism_1e20.txt"
nh=2.5 #2.5*4e19 = 1e20
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuation" -A "$nh"
times

# TEST 7b: With more ISM
out_file="$out_dir/jul6_nsxhnew_ism_4e20.txt"
nh=10 #10*4e19 = 4e20
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuation" -A "$nh"
times 

# TEST 8: With Effective Area applied
out_file="$out_dir/jul6_nsxhnew_ism_area.txt"
nh=5 #5*4e19 = 2e20
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuation" -A "$nh" -R "$inst_res"
times 

# TEST 9: With McPHAC Hydrogen
out_file="$out_dir/jul6_mcphacc_default.txt"
beaming=10
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 9a: With NSX Hydrogen Partially Ionized
out_file="$out_dir/jul6_nsxhpi_default.txt"
beaming=16
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times

# TEST 9b: With NSX Helium
out_file="$out_dir/jul6_nsxhe_default.txt"
beaming=15
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back"
times