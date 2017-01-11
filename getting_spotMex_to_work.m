clear, close all

cd /Users/kitung/Desktop/thesis_material/Spot
mex spotMex_trial.cpp -L/Users/kitung/Desktop/thesis_material/Spot -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx Instru.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.c
cd /Users/kitung/Desktop/thesis_material/Spot/ML2015_data/
obsdata = load('newa-10-err.txt');
cd /Users/kitung/Desktop/thesis_material/Spot


%setting up
auxOutput = cell(1,5);
%parameters

mass = 1.6;
radius = 11.8;
freq = 600;
inclination = 90;
emission = 90;
timeShift = 0;
numbins = 16;
modelchoice = 2;
rho = 0.4363;
spot_temperature = 2;
distance = 10;
numtheta = 4;
spectral_model = 2;
numbands = 30;
E_band_lower_1 = 3.5;
E_band_upper_1 = 12.5;
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
spotshape = 2;
obstime = 1.0405;
inst_curve = 0;
attenuation = 0;


%background = zeros(1,numbands);
background = 0.0705*ones(1,numbands);


cmd = '[Fspot,auxOutput{1}] = spotMex_trial(mass, radius, freq, inclination, emission, timeShift, numbins, modelchoice, rho, spot_temperature, distance, numtheta, spectral_model, numbands, E_band_lower_1, E_band_upper_1, beaming_model, spots_2, obsdata(:,1), bend_file_is, mr, b, psi, dcosa, toa, spotshape, obstime, inst_curve, attenuation';
for i = 1:numbands
    i
    cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
end
cmd = [cmd,');'];
disp(cmd)
eval(cmd);
                   

disp(Fspot)

