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
numtheta=1     # number of theta bins; for a small spot, only need one
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
rho=1e-2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.1 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters
E1Lower=0.3   # lower energy boundary
E1Uower=0.8   # upper energy boundary


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1 = small spot, BB
echo "start test 1"
out_file="$out_dir/bbsmall.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 1 = small spot, Graybody
echo "start test 2"
out_file="$out_dir/gsmall.txt"
beaming=1
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 3 = small spot, hydrogen
echo "start test 3"
out_file="$out_dir/hsmall.txt"
beaming=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 4 = small spot, helium
echo "start test 4"
out_file="$out_dir/hesmall.txt"
beaming=4
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

rho=1 # large spot;
numtheta=10 # number of theta bins;

# TEST 5 = large spot, BB
echo "start test 5"
out_file="$out_dir/bblarge.txt"
beaming=0
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 6 = large spot, Graybody
echo "start test 6"
out_file="$out_dir/glarge.txt"
beaming=1
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 7 = large spot, hydrogen
echo "start test 7"
out_file="$out_dir/hlarge.txt"
beaming=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 8 = large spot, helium
echo "start test 8"
out_file="$out_dir/helarge.txt"
beaming=4
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times