#!/bin/bash

# Scripts to run NICER code tests -- Sharon's computer settings
times
base="/Users/sharon/code/Albert"
exe_dir="$base/Spot-master-8"
#pwd
make spot
times
out_dir="$exe_dir/OBLATE-gamma"


normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=1     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numbands=1     # energy bands
spectraltype=0 # Blackbody observed at one Energy
beaming=0      # isotropic emission, no beaming
spotmodel=0


# doubles -- using "e" notation
spin=600 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1e-2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.35 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST OS1a = 600 Hz, tiny spot; oblate
numtheta=1 # number of theta bins; for a small spot, only need one
spotmodel=2
NS_model=2
out_file="$out_dir/test1a.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
times

# TEST OS1b = 600 Hz, big spot
NS_model=1
numbins=128
numtheta=100
rho=1.0
out_file="$out_dir/os1b.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1z = 600 Hz, big spot
NS_model=1
numbins=128
numtheta=1
rho=1e-2
inclination=90
emission=85
spin=600
out_file="$out_dir/os1z.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times


# TEST OS1c = 200 Hz, tiny spot; oblate
numtheta=1 # number of theta bins; for a small spot, only need one
rho=1e-2  # in radians
spin=200
out_file="$out_dir/os1c.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1d = 1 Hz, large spot; oblate
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=1.0
out_file="$out_dir/os1d.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST SP1d = 1 Hz, large spot; sphere
NS_model=3
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=1.0
out_file="$out_dir/sp1d.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1e = 600 Hz, large spot; oblate
numbins=128
NS_model=1
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=30 # in degrees
emission=60 # in degrees
out_file="$out_dir/os1e.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1f = 600 Hz, large spot; oblate
NS_model=1
numbins=128
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=80 # in degrees
emission=20 # in degrees
beaming=0
out_file="$out_dir/os1f.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1g = 600 Hz, large spot; oblate; cos^2 beaming
NS_model=1
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=30 # in degrees
emission=60 # in degrees
beaming=5
out_file="$out_dir/os1g.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1h = 600 Hz, large spot; oblate; sin^2 beaming
NS_model=1
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=30 # in degrees
emission=60 # in degrees
beaming=6
out_file="$out_dir/os1h.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1i = 600 Hz, large spot; oblate
NS_model=1
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=80 # in degrees
emission=20 # in degrees
beaming=5
out_file="$out_dir/os1i.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1j = 600 Hz, large spot; oblate
NS_model=1
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=80 # in degrees
emission=20 # in degrees
beaming=6
out_file="$out_dir/os1j.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1e = 600 Hz, large spot; spherical; circular spot in rotating frame
numbins=128
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=30 # in degrees
emission=60 # in degrees
spotmodel=1
beaming=0
out_file="$out_dir/sp1e-circ.txt"
NS_model=3
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1e = 600 Hz, large spot; oblate; circular spot in rotating frame
NS_model=1
out_file="$out_dir/os1e-circ.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1f = 600 Hz, large spot; oblate
NS_model=1
numbins=128
numtheta=40 # number of theta bins; for a small spot, only need one
rho=1  # in radians
spin=600
inclination=80 # in degrees
emission=20 # in degrees
beaming=0
spotmodel=1
out_file="$out_dir/os1f-circ.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OSALL = OBLATE STAR - Spot covers all of star.
#NS_model=1
numbins=512
numtheta=100 # number of theta bins; for a small spot, only need one
rho=3.1415926535  # in radians
spin=1
inclination=90 # in degrees
emission=0 # in degrees
beaming=0
spotmodel=0
NS_model=1
spin=1
out_file="$out_dir/osall-1.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times
NS_model=1
spin=200
out_file="$out_dir/osall-200.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times
NS_model=1
spin=400
out_file="$out_dir/osall-400.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times
NS_model=1
spin=600
out_file="$out_dir/osall-600.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# integers
numtheta=1     # number of theta bins; for a small spot, only need one
NS_model=1     # 1 (oblate) or 3 (spherical); don't use 2
numbins=128    # phase bins
numbands=1     # energy bands
beaming=0      # isotropic emission, no beaming
spotmodel=0

# doubles -- using "e" notation
spin=600 # in Hz
mass=1.4 # in Msun
radius=12 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1e-2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=0.35 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters

# TEST OS1a = 600 Hz, tiny spot; oblate
spectraltype=0 # Blackbody observed at one Energy
out_file="$out_dir/os1a.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel"
#times

# TEST OS1k = 600 Hz, tiny spot; oblate
spectraltype=1 # Funny NICER Line flux
numbands=300
rho=1e-2
spin=400 # in Hz
l1=0.995 # in keV
l2=1.005 # in keV
e0=0.7   # in keV
de=1e-3  # in keV
out_file="$out_dir/os1k.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel" -v "$l1" -V "$l2" -x "$e0" -X "$de"
#times

# TEST OS1k = 600 Hz, tiny spot; oblate
spectraltype=1 # Funny NICER Line flux
numbands=300
rho=1.0
numtheta=100
spin=400 # in Hz
l1=0.995 # in keV
l2=1.005 # in keV
e0=0.7   # in keV
de=1e-3  # in keV
out_file="$out_dir/os1l.txt"
## RUNNING THE CODE
#./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$beaming" -s "$spectraltype" -S "$numbands" -P "$spotmodel" -v "$l1" -V "$l2" -x "$e0" -X "$de"
#times