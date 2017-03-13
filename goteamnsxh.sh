#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/kitung/Desktop/thesis_material"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/team_hw"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=4     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=32    # phase bins
numbands=15     # energy bands
spectraltype=2 # BB Integrated energies
elo=0.1
ehi=3.1
beaming=0      # bb
spotmodel=0    # circular in the static frame, no gamma
obstime=1000000.0      #length of observation (in seconds)
back=0.0    #constant background (same in all energy bands)

# doubles -- using "e" notation
spin=0.01 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.1  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.09669 # in keV (Surface temperature)
distance=0.2 # distance in kpc

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1: NSX Hydrogen, 1.1e6K, 90/90, monochromatic
out_file="$out_dir/nsxh_6.05_9090_mono.txt"
beaming=5      # NSX hydrogen
temp=0.09669 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 2: NSX Hydrogen, 1.1e6K, 90/90, integrated
out_file="$out_dir/nsxh_6.05_9090_inte.txt"
beaming=5      # NSX hydrogen
temp=0.09669 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times

# TEST 3: NSX Hydrogen, 1.1e6K, 30/80, monochromatic
out_file="$out_dir/nsxh_6.05_3080_mono.txt"
beaming=5      # NSX hydrogen
temp=0.09669 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 4: NSX Hydrogen, 1.1e6K, 30/80, integrated
out_file="$out_dir/nsxh_6.05_3080_inte.txt"
beaming=5      # NSX hydrogen
temp=0.09669 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times

# TEST 5: NSX Hydrogen, 2e6K, 90/90, monochromatic
out_file="$out_dir/nsxh_6.3_9090_mono.txt"
beaming=5      # NSX hydrogen
temp=0.1758 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 6: NSX Hydrogen, 2e6K, 90/90, integrated
out_file="$out_dir/nsxh_6.3_9090_inte.txt"
beaming=5      # NSX hydrogen
temp=0.1758 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times

# TEST 7: NSX Hydrogen, 2e6K, 30/80, monochromatic
out_file="$out_dir/nsxh_6.3_3080_mono.txt"
beaming=5      # NSX hydrogen
temp=0.1758 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 8: NSX Hydrogen, 2e6K, 30/80, integrated
out_file="$out_dir/nsxh_6.3_3080_inte.txt"
beaming=5      # NSX hydrogen
temp=0.1758 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times

# TEST 9: NSX Hydrogen, 5e5K, 90/90, monochromatic
out_file="$out_dir/nsxh_5.69_9090_mono.txt"
beaming=5      # NSX hydrogen
temp=0.042 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 10: NSX Hydrogen, 5e5K, 90/90, integrated
out_file="$out_dir/nsxh_5.69_9090_inte.txt"
beaming=5      # NSX hydrogen
temp=0.042 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=90 # in degrees
emission=90 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times

# TEST 11: NSX Hydrogen, 5e5K, 30/80, monochromatic
out_file="$out_dir/nsxh_5.69_3080_mono.txt"
beaming=5      # NSX hydrogen
temp=0.042 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=0 # Monochromatic
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res" 
times

# TEST 12: NSX Hydrogen, 5e5K, 30/80, integrated
out_file="$out_dir/nsxh_5.69_3080_inte.txt"
beaming=5      # NSX hydrogen
temp=0.042 # in keV (Surface temperature)
inst_res=2     # 1 cm^2 for all bands
inclination=30 # in degrees
emission=80 # in degrees
spectraltype=3 # Integrated
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -R "$inst_res"
times