#!/bin/bash
# January 19, 2018 - Version to test 10^6+10^6 Synthetic Data
# Scripts to run NICER code tests -- Sharon's computer settings
times

base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-29"

make spot
times
out_dir="$exe_dir/outmr1e6"

# integers
NS_model=1     # 1 (oblate) or 3 (spherical); 2 is for old-fashioned shape model
numbins=16     # phase bins
numbands=301   # energy bands
spectraltype=0 # 
beaming=10     # McPhac
spotmodel=0    # circular in the static frame, no gamma
inst_res=0
attenuation=0  # 

# doubles -- using "e" notation
spin=600 # in Hz
mass=1.4 # in Msun
radius=12.0 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1.0  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.231139 # in keV (Surface temperature)
distance=0.2 # distance in kpc
elo=0.095
ehi=3.105     
#obstime=1.0 #length of observation (in seconds)
#obstime=41.1756 #length of observation (in seconds) Suitable for 10^6 counts
obstime=41.176 #length of observation (in seconds) Suitable for 10^6 counts


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi


# TEST 2:  Instrumental Response added
inst_res=1
spectraltype=0
beaming=10
numtheta=40
out_file="$out_dir/spot6.txt"
#data_file="Slavko/CU_high_accuracy_1spot_10k_plaw_10k_poisson_sampled.txt"
data_file="Slavko/new1e6.txt"

radius=12.0


## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "output.txt" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_file" -R "$inst_res" -I "$data_file" -K "Background/background6a_guess.txt"
times


