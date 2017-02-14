function extPar = init  %defining an external parameter sent to the executable; equivalent of command line arguments to main

%spotDir=fileparts(which(mfilename));
spotDir=fileparts(which('init.m'));
%
% This version of init.m is for running a comparison with Slavko's
% synthetic data
% ------------
% Compile the functions...
%
% If you see bizarre errors when mex evaluates..
% There is a file called mexopts.sh that needs to be modified.
% It is located in /home/<your user name>/.matlab/R2010a/mexopts.sh (may depend on MATLAB version)
% look in there and you will see about four places where "-ansi" is specified as a flag.
% Simply remove the "-ansi" and it should compile.
%
cmd=['mex spotMex_trial.cpp -L', spotDir, ' -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp'];
cd(spotDir);
disp('------ Init ------');
disp(cmd);
which('spotMex_trial.cpp');
eval(cmd);
%
%
% ------------
% extPar is external parameters
%
% Parameters for comparison with ML2015
%
%extPar.fixed.mass=1.6;                % mass in MSun
%extPar.fixed.radius=11.8;               % radius in km
extPar.fixed.freq=300;                % spin frequency of NS, in Hz
%extPar.fixed.inclination=90;          % spot inclination angle in degrees
%extPar.fixed.emission=90;             % spot emission angle in degrees
%extPar.fixed.phaseshift=0;            % spot phaseshift angle in degrees
extPar.fixed.numbins=32;
extPar.fixed.modelchoice=1;           % 1 (oblate 2014 A&M), 2 (oblate 2007 MLCB), or 3 (spherical);
extPar.fixed.rho=0.17453;                   % angular radius of the emitting hot spot, in radians
extPar.fixed.spot_temperature=0.35;      % temperature of the hot spot in frame of the NS, in keV
extPar.fixed.distance=0.3;              % distance from us to the NS, in kpc
extPar.fixed.numtheta=10;              % number of theta bins.
extPar.fixed.spectral_model=2;		% 2 is for integrated bb and variation, 3 is for atmosphere flux integrated within energy bands
extPar.fixed.numbands=15;              % 
extPar.fixed.E_band_lower_1=0.1;      % lower bound of first energy band, in keV
extPar.fixed.E_band_upper_1=3.1;      % upper bound of first energy band, in keV
extPar.fixed.beaming_model=7;         % 0 for bb, 1 for bb+chandra gray, 2 for bb+hopf gray, 3 for hydrogen, 4 for helium
extPar.fixed.spots_2=1;				% 1 for 1 spot, 2 for 2 spots
extPar.fixed.obstime=0.2;        % Time (in seconds) that source was observed
extPar.fixed.inst_curve=0;
extPar.fixed.attenuation=0;
extPar.fixed.bend_file_is=1;          % 1 means bend file is read
if extPar.fixed.bend_file_is == 1
    bendfile = load('angles1000-1.txt');
    bendfile1 = reshape(bendfile(:,6),301,1001);
    extPar.fixed.mr = bendfile1(1,:)';
    extPar.fixed.b = bendfile(:,2);
    extPar.fixed.psi = bendfile(:,3);
    extPar.fixed.dcosa = bendfile(:,4);
    extPar.fixed.toa = bendfile(:,5);
else
    extPar.fixed.mr = 0;
    extPar.fixed.b = 0;
    extPar.fixed.psi = 0;
    extPar.fixed.dcosa = 0;
    extPar.fixed.toa = 0;
end
extPar.fixed.spotshape=0;
%
%
% ------------
% Load data.
% NEED TO CHANGE THIS IN outputFerret.m TOO, ~ line 80
exe_dir='/home/kitung/Spot';
% fileName=strcat(exe_dir,'/input/bb_only_fake_data.txt'); % *** need to change this in outputFerret too
fileName=strcat(exe_dir,'/input/sbpoisson3.txt'); % *** need to change this in outputFerret too
eval('dataFile=load(fileName);');
% ------------
% Process data
% 
extPar.obsdata2.t=dataFile(:,1)';
for i = 1:extPar.fixed.numbands
    extPar.obsdata2.f(i,:)=dataFile(:,2*i)';
    extPar.obsdata2.err(i,:)=dataFile(:,2*i+1)';
end
extPar.fixed.background = zeros(1,extPar.fixed.numbands);
%extPar.obsdata2.numbins=size(dataFile,1);
end
