%clear, close all

tic;
disp('Compiling!')


cd /home/kitung/Spot

mex spotMex_new.cpp -L/home/kitung/Spot -cxx Units.cpp -cxx Chi.cpp -cxx McPhac.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp -cxx TimeDelays.cpp -cxx BlackBody.cpp

time=toc
disp('Setting Things up!')
tic;

%Load in the Synthetic Data
disp('Loading the Synthetic Data')
%obsdata = load('Slavko/CU_high_accuracy_1spot_10k_plaw_10k_poisson_sampled.txt');
obsdata = load('synthetic_2e6_realisation.dat');
datatime = obsdata(:,2);
dataflux = obsdata(:,3);
obsdatatime = reshape(datatime,32,300);
obsdatanew = reshape(dataflux,32,300);

%Load the Instrumental Response Matrix
disp('Loading the Instrumental Response Matrix')
instru = load('Area/NICER_May2014_rsp.txt');

%Load the Hydrogen Atmosphere 
disp('Loading the Hydrogen Atmosphere Model')
atmotable=fopen('atmosphere/Hatm8000dT0.05.bin');
a = fread(atmotable,'double');   
fclose(atmotable);
atmo = a([1:1595000]*5);
angl = a([1:1595000]*5-1);
energytable = load('atmosphere/Energy.txt');
energy = energytable(:,2);

auxOutput = cell(1,5);
%parameters

disp('Setting the Parameter Values')
mass = 1.4;
radius = 12.0;
freq = 600;
inclination = 90;
emission = 90;
timeShift = 0;
numbins = 16;
modelchoice = 1; %oblate
rho = 1.0;
spot_temperature = 0.231139;
distance = 0.2;
numtheta = 6;
spectral_model = 0;
numbands = 301;
E_band_lower_1 = 0.095;
E_band_upper_1 = 3.105;
beaming_model = 10;
spots_2 = 1;
bendfile = load('angles1000-1.txt');
bend_file_is = 1;
bendfile1 = reshape(bendfile(:,6),301,1001);
mr = bendfile1(1,:)';
b = bendfile(:,2);
psi = bendfile(:,3);
dcosa = bendfile(:,4);
toa = bendfile(:,5);
spotshape = 0;
%obstime = 41.176;
obstime = 28.56;
inst_curve = 1;
attenuation = 5;
nh = 1.1;  % in units of 10^{18} cm^2
inte = atmo;

disp('Loading the guess for the background')
%background = zeros(1,numbands);
background = load('Background/background6a.txt');

disp('Loading the ISM Absorption Model')
% 3 column file format: col 1=NH x 10^{18} cm^2, col 2=Energy,
% col 3=attenuation factor
atten_full = load('ISM/tbnew_full.txt');
atten = atten_full(:,3);

time=toc

cmd = '[Fspot,auxOutput{1}] = spotMex_new(mass, radius, freq, inclination, emission, timeShift, numbins, modelchoice, rho, spot_temperature, distance, numtheta, spectral_model, numbands, E_band_lower_1, E_band_upper_1, beaming_model, spots_2, obsdatatime(:,1), bend_file_is, mr, b, psi, dcosa, toa, spotshape, obstime, inst_curve, attenuation, inte, angl, energy';
for i = 1:numbands-1
    
    %cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
    cmd = [cmd,', obsdatanew(:,',num2str(i),'),background(',num2str(i),')'];
    
end

for i = 1:400
    cmd = [cmd,', instru(',num2str(i),',:)'];
end

cmd = [cmd,',atten, nh'];

cmd = [cmd,');'];
disp(cmd)

tic;
eval(cmd);
time=toc                   

disp(Fspot)

