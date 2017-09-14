#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
#base="/Users/kitung/Desktop/thesis_material"
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-26"
#pwd
make spot
times
out_dir="$exe_dir/Sept8"

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
obstime=1.0 #length of observation (in seconds)

## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi





# TEST 2:  Instrumental Response added
inst_res=1
spectraltype=0
beaming=10
numtheta=6
out_file="$out_dir/intspot.txt"

## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -E "$deltatheta" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "output.txt" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -u "$elo" -U "$ehi" -P "$spotmodel" -b "angles1000.txt" -Z "$obstime"  -O "$out_file" -R "$inst_res"
times


cd "$out_dir"



    slavko="CU_highAccuracy_1spot.txt"
    approx="intspot.txt"
    poisson="Slavko/CU_highAccuracy_1spot_pl_poisson_sampled_01.txt"



    DiffLogLikely -i "$poisson" -j "$approx" -k "$slavko" -o "slavdiff$n.txt"




#DiffCurve -i "CU_highAccuracy_1spot.txt" -j "$out_file" -o "intdiff_$numtheta.txt"

#DiffCurve -i "CU_high_accuracy_1spot_10k_plaw_10k_poisson_sampled.txt" -j "$out_file" -o "lowdiff_$numtheta.txt"


#DiffCurve -i "spotnofoldunsamp_20170725-im.txt" -j "$out_file" -o "IMdiff160-256.txt"
#DiffCurve -i "spotnofoldunsamp_20170725-im.txt" -j "$out_file" -o "AMSdiff40.txt"
#DiffCurve -i "CU_1spot_unsampled_noresponse_20170726.txt" -j "$out_file" -o "diff_noresp.txt"

#DiffCurve -i "bbspotnofoldunsamp-ab1.txt" -j "$out_file" -o "bbdiff80.txt"

#DiffCurve -i "spotnofold-ab6.txt" -j "$out_file" -o "olddiff.txt"