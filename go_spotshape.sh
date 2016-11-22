#!/bin/bash

# Scripts to run NICER code tests
times
base="/Users/kitung/Desktop/thesis_material/"
exe_dir="$base/Spot"
#pwd
make spot
times
out_dir="$exe_dir/SpotShape"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=3     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
#numphi=20	   # phi bins	
numbands=2     # energy bands
spectraltype=3 # 
beaming=3      # Helium


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

#input fake observation
fakedata="he.txt"

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1 = one small spot, does not go over pole
echo "start test 1"
out_file="$out_dir/1_sml_no.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming"
times

# TEST 2 = two small spots, does not go over pole
echo "start test 2"
out_file="$out_dir/2_sml_no.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 3 = one large spot, does not go over pole
echo "start test 3"
rho=1
numtheta=5     # number of theta bins; for a small spot, only need one
out_file="$out_dir/1_lrg_no.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming"
times

# TEST 4 = two large spots, does not go over pole
echo "start test 4"
out_file="$out_dir/2_lrg_no.txt"
numtheta=
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times





# TEST 5 = small spot, goes over pole
echo "start test 5"
rho=1e-2
numtheta=1     # number of theta bins; for a small spot, only need one
inclination=90
emission=0.5
out_file="$out_dir/sml_pole.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times

# TEST 6 = large spot, goes over pole
echo "start test 6"
emission=10
rho=1
numtheta=5     # number of theta bins; for a small spot, only need one
out_file="$out_dir/lrg_pole.txt"
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -s "$spectraltype" -S "$numbands" -u "$E1Lower" -U "$E1Upper" -g "$beaming" -2
times
