clear, close all

tic;
disp('Compiling!')

%cd /Users/kitung/Desktop/thesis_material/Spot
cd /home/kitung/Spot
mex spotMex_trial.cpp -L/home/kitung/Spot -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp
cd /home/kitung/Spot/input/

obsdata = load('smmpoisson1.txt');
%cd /Users/kitung/Desktop/thesis_material/Spot
cd /home/kitung/Spot

time=toc

tic;
disp('Setting Things up!')

%setting up
auxOutput = cell(1,5);
%parameters

mass = 1.6;
radius = 11.8;
freq = 300;
inclination = 90;
emission = 90;
timeShift = 0;
numbins = 32;
modelchoice = 1;
rho = 0.17453;
spot_temperature = 0.35;
distance = 0.3;
numtheta = 18;
spectral_model = 2;
numbands = 15;
E_band_lower_1 = 0.1;
E_band_upper_1 = 3.1;
beaming_model = 7;
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
obstime = 6.4;
inst_curve = 0;
attenuation = 0;


%background = zeros(1,numbands);
background = 0.005*ones(1,numbands);

time=toc

cmd = '[Fspot,auxOutput{1}] = spotMex_trial(mass, radius, freq, inclination, emission, timeShift, numbins, modelchoice, rho, spot_temperature, distance, numtheta, spectral_model, numbands, E_band_lower_1, E_band_upper_1, beaming_model, spots_2, obsdata(:,1), bend_file_is, mr, b, psi, dcosa, toa, spotshape, obstime, inst_curve, attenuation';
for i = 1:numbands
    
    cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
end
cmd = [cmd,');'];
disp(cmd)
tic;
eval(cmd);
time=toc                   

disp(Fspot)

