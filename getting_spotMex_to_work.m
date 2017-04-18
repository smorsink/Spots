clear, close all

tic;
disp('Compiling!')

cd /Users/kitung/Desktop/thesis_material/Spot
%cd /home/kitung/Spot
mex spotMex_trial.cpp -L/Users/kitung/Desktop/thesis_material/Spot -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.cpp

time=toc
disp('Setting Things up!')
tic;

cd /Users/kitung/Desktop/thesis_material/Spot/input/
%cd /home/kitung/Spot/input/
obsdata = load('mcphac_apr_11_uncert.txt');
%obsdata = load('mcphac_apr_12_withnoise.txt');
%cd /Users/kitung/Desktop/thesis_material/Spot/atmosphere/mcphacc
cd /Users/kitung/Desktop/thesis_material/atmospheres/cole_mcphac
atmotable=fopen('Hatm8000dT0.05.bin');
a = fread(atmotable,'double');   
fclose(atmotable);
atmo = a([1:1595000]*5);
%atmo = load('cole_mcphac_table2.txt');
cd /Users/kitung/Desktop/thesis_material/Spot
%cd /home/kitung/Spot

%setting up
auxOutput = cell(1,5);
%parameters

mass = 1.48;
radius = 11.8;
freq = 300;
inclination = 90;
emission = 90;
timeShift = 0;
numbins = 32;
modelchoice = 3;
rho = 0.2;
spot_temperature = 0.09669;
distance = 0.3;
numtheta = 18;
spectral_model = 3;
numbands = 15;
E_band_lower_1 = 0.1;
E_band_upper_1 = 3.1;
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
obstime = 1;
inst_curve = 0;
attenuation = 0;
%inte = atmo(:,3);
inte = atmo;


%background = zeros(1,numbands);
background = 0.001*ones(1,numbands);

time=toc

cmd = '[Fspot,auxOutput{1}] = spotMex_trial(mass, radius, freq, inclination, emission, timeShift, numbins, modelchoice, rho, spot_temperature, distance, numtheta, spectral_model, numbands, E_band_lower_1, E_band_upper_1, beaming_model, spots_2, obsdata(:,1), bend_file_is, mr, b, psi, dcosa, toa, spotshape, obstime, inst_curve, attenuation, inte';
for i = 1:numbands
    
    cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
end
cmd = [cmd,');'];
disp(cmd)
tic;
eval(cmd);
time=toc                   

disp(Fspot)

