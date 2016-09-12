#!/bin/bash

# Scripts to run NICER code tests
times
base="/Users/kitung/Desktop/thesis_material/"
exe_dir="$base/Spot-master"
#pwd
make spot
times
out_dir="$exe_dir/ebands"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=10     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numphi=20	   # phi bins	
numbands=5     # energy bands
spectraltype=2 # Integrated flux with energy bands
beaming=0      # Blackbody, no beaming


# doubles
spin=400 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.1 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters
E1Lower=0.3   # lower energy boundary
E1Upper=0.8   # upper energy boundary


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi


# TEST 1 = BB, observed from 0.3-0.8 keV
echo "start test 1"
out_file="$out_dir/bblow.txt"
beaming=0
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 2 = Graybody, observed 0.3-0.8 keV
echo "start test 2"
out_file="$out_dir/glow.txt"
beaming=1
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 3 = hydrogen, observed 0.3-0.8 keV
echo "start test 3"
out_file="$out_dir/hlow.txt"
beaming=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 4 = helium, observed 0.3-0.8 keV
echo "start test 4"
out_file="$out_dir/helow.txt"
beaming=4
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

E1Lower=1   # lower energy boundary
E1Upper=1.5   # upper energy boundary

# TEST 1 = BB, observed from 1-1.5 keV
echo "start test 5"
out_file="$out_dir/bbhigh.txt"
beaming=0
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 2 = Graybody, observed from 1-1.5 keV
echo "start test 6"
out_file="$out_dir/ghigh.txt"
beaming=1
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 3 = hydrogen, observed from 1-1.5 keV
echo "start test 7"
out_file="$out_dir/hhigh.txt"
beaming=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 4 = helium, observed from 1-1.5 keV
echo "start test 8"
out_file="$out_dir/hehigh.txt"
beaming=4
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times


