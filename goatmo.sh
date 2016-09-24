#!/bin/bash

# Scripts to run NICER code tests
times
base="/Users/kitung/Desktop/thesis_material/"
exe_dir="$base/Spot-master"
#pwd
make spot
times
out_dir="$exe_dir/ATMOTEST"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numphi=20	   # phi bins	
numbands=5     # energy bands
spectraltype=2 # Integrated flux with energy bands
#spectraltype=3 # nothing calculated; time is to load atmosphere only.
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
E1Upper=0.8   # upper energy boundary


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi


# TEST 1 = hydrogen, observed 0.3-0.8 keV
echo "start test 1"
out_file="$out_dir/h.txt"
beaming=3
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 2 = helium, observed 0.3-0.8 keV
echo "start test 2"
out_file="$out_dir/he.txt"
beaming=4
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 3 = helium2, observed 0.3-0.8 keV
# 6, 5, 4, 3, 3 steps for each of the five bands, respectively
echo "start test 3"
out_file="$out_dir/he2.txt"
spectraltype=3
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 4"
# 0 step
out_file="$out_dir/he2-1.txt"
numbands=1
E1Lower=0.7   # lower energy boundary
E1Upper=0.71   # upper energy boundary
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 5"
# 1 step
out_file="$out_dir/he2-2.txt"
E1Lower=0.7   # lower energy boundary
E1Upper=0.73   # upper energy boundary
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 5"
# 2 step
out_file="$out_dir/he2-3.txt"
E1Lower=0.7   # lower energy boundary
E1Upper=0.76   # upper energy boundary
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 6"
# 0 step
out_file="$out_dir/he-1.txt"
spectraltype=2
E1Lower=0.7   # lower energy boundary
E1Upper=0.71   # upper energy boundary
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 7"
# 1 step
out_file="$out_dir/he-2.txt"
E1Lower=0.7   # lower energy boundary
E1Upper=0.73   # upper energy boundary
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

echo "start test 8"
# 2 step
out_file="$out_dir/he-3.txt"
E1Lower=0.7   # lower energy boundary
E1Upper=0.76   # upper energy boundary
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -P "$numphi" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times
