function extPar = init  %defining an external parameter sent to the executable; equivalent of command line arguments to main

% This version of the code reads in a fixed background. Only geometric
% parameters are allowed to vary

spotDir=fileparts(which('init.m'));
%extPar.dataDir='/home/kitung/Ferret_Runs/May12a';
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
%cmd=['mex spotMex_trial.cpp -L', spotDir, ' -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp'];

disp('----  Hello ------');
cd(spotDir);
mex spotMex_new.cpp -L/home/kitung/Spot -cxx Units.cpp -cxx Chi.cpp -cxx McPhac.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp -cxx TimeDelays.cpp -cxx BlackBody.cpp

%cd(spotDir);
disp('------ Init ------');

%
% ------------
% extPar is external parameters
%
% Parameters for comparison with smmpoisson1.txt
%
extPar.fixed.mass=1.4;                % mass in MSun
extPar.fixed.radius=12.0;               % radius in km
extPar.fixed.freq=600;                % spin frequency of NS, in Hz
extPar.fixed.inclination=60;          % spot inclination angle in degrees
extPar.fixed.emission=100;             % spot emission angle in degrees
extPar.fixed.phaseshift=0;            % spot phaseshift angle in degrees
extPar.fixed.numbins=32;
extPar.fixed.modelchoice=1;           % 1 (oblate 2014 A&M), 2 (oblate 2007 MLCB), or 3 (spherical);
extPar.fixed.rho=1.0;                   % angular radius of the emitting hot spot, in radians
extPar.fixed.spot_temperature=0.231139;      % temperature of the hot spot in frame of the NS, in keV
extPar.fixed.distance=0.2;              % distance from us to the NS, in kpc
extPar.fixed.numtheta=40;              % number of theta bins.
extPar.fixed.spectral_model=0;		% 2 is for integrated bb and variation, 3 is for atmosphere flux integrated within energy bands
extPar.fixed.numbands=301;              % 
extPar.fixed.E_band_lower_1=0.095;      % lower bound of first energy band, in keV
extPar.fixed.E_band_upper_1=3.105;      % upper bound of first energy band, in keV
extPar.fixed.beaming_model=10;         % 0 for bb, 1 for bb+chandra gray, 2 for bb+hopf gray, 3 for hydrogen, 4 for helium
extPar.fixed.spots_2=1;				% 1 for 1 spot, 2 for 2 spots
%extPar.fixed.obstime=41.176;        % Time (in million seconds) that source was observed
%extPar.fixed.obstime=28.56;
extPar.fixed.obstime=24.913;
extPar.fixed.inst_curve=1;
extPar.fixed.attenuation=5;
extPar.fixed.nh=1;
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

disp(extPar.fixed.spotshape);

extPar.fixed.instru = load('Area/NICER_May2014_rsp.txt');
disp('Response Maxtrx Loaded');

%Load Atmosphere
%cd /home/kitung/atmospheres/cole_mcphac
atmotable=fopen('atmosphere/Hatm8000dT0.05.bin');
a = fread(atmotable,'double');   
fclose(atmotable);
atmo = a([1:1595000]*5);
angl = a([1:1595000]*5-1);
extPar.fixed.inte = atmo;
extPar.fixed.angl = angl(1:50);
energytable = load('atmosphere/Energy.txt');
energy = energytable(:,2);
extPar.fixed.energy=energy;
%disp(energy);
%
disp('Loading the ISM Absorption Model')
% 3 column file format: col 1=NH x 10^{18} cm^2, col 2=Energy,
% col 3=attenuation factor
atten_full = load('ISM/tbnew_full.txt');
extPar.fixed.atten = atten_full(:,3);
%
% ------------
% Load data.
%obsdata = load('Slavko/CU_high_accuracy_1spot_10k_plaw_10k_poisson_sampled.txt');
%obsdata = load('Slavko/new1e6.txt');
obsdata = load('synthetic_2e6_realisation.dat');
datatime = obsdata(:,2);
dataflux = obsdata(:,3);

obsdatatime = reshape(datatime,32,300);
obsdatanew = reshape(dataflux,32,300);


extPar.fixed.obsdata.t=obsdatatime;
for i = 1:extPar.fixed.numbands-1
    extPar.fixed.obsdata.f(:,i)=obsdatanew(:,i);
 %   extPar.obsdata2.f(i,:)=dataFile(:,2*i)';
 %   extPar.obsdata2.err(i,:)=dataFile(:,2*i+1)';
end

% Load Background
%extPar.fixed.background = 0.0*ones(1,extPar.fixed.numbands);
%extPar.fixed.background = zeros(1,extPar.fixed.numbands);
extPar.fixed.background = load('Background/background60-110.txt');
%extPar.obsdata2.numbins=size(dataFile,1);

%disp(extPar.fixed.mass)

cmd =  '[Fspot,auxOutput{1}]=spotMex_new(extPar.fixed.mass, extPar.fixed.radius, extPar.fixed.freq, extPar.fixed.inclination, extPar.fixed.emission, extPar.fixed.phaseshift, extPar.fixed.numbins, extPar.fixed.modelchoice, extPar.fixed.rho, extPar.fixed.spot_temperature, extPar.fixed.distance, extPar.fixed.numtheta, extPar.fixed.spectral_model, extPar.fixed.numbands, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.fixed.spots_2, extPar.fixed.obsdata.t, extPar.fixed.bend_file_is, extPar.fixed.mr, extPar.fixed.b, extPar.fixed.psi, extPar.fixed.dcosa, extPar.fixed.toa, extPar.fixed.spotshape, extPar.fixed.obstime, extPar.fixed.inst_curve, extPar.fixed.attenuation, extPar.fixed.inte, extPar.fixed.angl, extPar.fixed.energy';
for j = 1:300    
    %cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
    cmd = [cmd,', extPar.fixed.obsdata.f(:,',num2str(j),'),extPar.fixed.background(',num2str(j),')'];
    % cmd = [cmd,', extPar.fixed.obsdata.f(:,',num2str(j),'),par.back(',num2str(j),')'];
end

for j = 1:400
    cmd = [cmd,', extPar.fixed.instru(',num2str(j),',:)'];
end
cmd = [cmd,',extPar.fixed.atten, extPar.fixed.nh'];
cmd = [cmd,');'];
tic;
        disp(cmd)
        disp(i);
        %disp(timeShift);
        eval(cmd);
        disp(Fspot);
        time=toc


end
