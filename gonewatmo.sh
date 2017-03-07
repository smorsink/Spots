#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
#base="/Users/kitung/Desktop/thesis_material"
base="/home/kitung"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/new_atmo_test"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=32    # phase bins
numbands=30     # energy bands
spectraltype=2 # BB Integrated energies
elo=0.05
ehi=3.05
beaming=0      # bb
spotmodel=0    # circular in the static frame, no gamma
obstime=1040500.0      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)

# doubles -- using "e" notation
spin=600 # in Hz
mass=1.48 # in Msun
radius=11.8 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.4363  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.09669 # in keV (Surface temperature)
#distance=10.0 # distance in kpc
#distance=5.833 # PSR J0030, 1.8e20 m in kpc
distance=1.296 # PSR J0437, 4e19 m in kpc

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1: BB
out_file="$out_dir/first_test_bb.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times

# TEST 2: NSATMOS
out_file="$out_dir/first_test_nsatmos.txt"
beaming=3      # hydrogen
spectraltype=3 # Atmosphere Integrated energies
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times

# TEST 3: NSX Helium
out_file="$out_dir/first_test_nsxhe.txt"
beaming=4      # hydrogen
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times

# TEST 4: NSATMOS, monochromatic
out_file="$out_dir/nsatmos_mono.txt"
beaming=3      # hydrogen
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times

# TEST 5: NSX Hydrogen, monochromatic
out_file="$out_dir/nsxh_mono.txt"
beaming=5      # NSX hydrogen
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times

# TEST 6: NSATMOS, integrated
out_file="$out_dir/nsatmos_inte.txt"
beaming=3      # hydrogen
spectraltype=3 # integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times


# TEST 7: NSX Hydrogen, integrated
out_file="$out_dir/nsxh_inte.txt"
beaming=5      # hydrogen
spectraltype=3 # integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -2
times
