#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/kitung/Desktop/thesis_material"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/instru_test"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=4     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=16    # phase bins
numbands=30     # energy bands
spectraltype=2 # Blackbody energy bands
elo=0.05
ehi=3.05
beaming=7      # Hopf function limb-darkening
spotmodel=2    # circular in the static frame, no gamma
obstime=1040500.0      #length of observation (in seconds)
back=0.0705    #constant background (same in all energy bands)

# doubles -- using "e" notation
spin=600 # in Hz
mass=1.6 # in Msun
radius=11.8 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=0.4363  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=2.0 # in keV (Surface temperature)
#distance=10.0 # distance in kpc
distance=5.833 # PSR J0030, 1.8e20 m in kpc
#distance=1.296 # PSR J0437, 4e19 m in kpc

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi


# TEST 1: no attenuation
out_file="$out_dir/J0030_none.txt"
attenuate=0
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuate"
#times


# TEST 2: attenuation by WABS, J0030
out_file="$out_dir/J0030_wabs.txt"
attenuate=1
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuate"
#times

# TEST 2: attenuation by TBABS, J0030
out_file="$out_dir/J0030_tbabs.txt"
attenuate=2
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuate"
#times


# TEST 3: background as listed by each energy band, J0030
out_file="$out_dir/J0030_backlist.txt"
back_file="J0030_fake_backlist"
attenuate=0
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -K "$back_file" -a "$attenuate"
#times

# TEST 4: background list having wrong number of bands
out_file="$out_dir/J0030_wronglist.txt"
back_file="J0030_wrong_backlist"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -K "$back_file" -a "$attenuate"
#times


# TEST 5: see if instrument response flag is working 

out_file="$out_dir/J0030_inst.txt"
attenuate=0
inst_res=1
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime" -k "$back" -a "$attenuate" -R "$inst_res"
#times


