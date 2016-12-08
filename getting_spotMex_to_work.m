clear, close all

cd /Users/kitung/Desktop/thesis_material/Spot
mex spotMex_trial.cpp -L/Users/kitung/Desktop/thesis_material/Spot -cxx Units.cpp -cxx Chi.cpp -cxx Atmo.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp -cxx nrutil.c -cxx interp.c
cd /Users/kitung/Desktop/thesis_material/Spot/input/
obsdata = load('bb_only_fake_data.txt');
cd /Users/kitung/Desktop/thesis_material/Spot


%setting up
auxOutput = cell(1,5);
%parameters

mass = 1.4;
radius = 12;
freq = 600;
inclination = 90;
emission = 90;
timeShift = 0;
numbins = 16;
modelchoice = 2;
rho = 1;
spot_temperature = 2;
distance = 3.08567758e20;
numtheta = 5;
spectral_model = 2;
numbands = 5;
E_band_lower_1 = 0.3;
E_band_upper_1 = 1.8;
beaming_model = 0;
spots_2 = 1;
bendfile = load('angles100-1.txt');
bend_file_is = 1;
bendfile1 = reshape(bendfile(:,6),151,101);
mr = bendfile1(1,:)';
b = bendfile(:,2);
psi = bendfile(:,3);
dcosa = bendfile(:,4);
toa = bendfile(:,5);
spotshape = 0;

%background = zeros(1,numbands);
background = 0.01*ones(1,numbands);


cmd = '[Fspot,auxOutput{1}] = spotMex_trial(mass, radius, freq, inclination, emission, timeShift, numbins, modelchoice, rho, spot_temperature, distance, numtheta, spectral_model, numbands, E_band_lower_1, E_band_upper_1, beaming_model, spots_2, obsdata(:,1), bend_file_is, mr, b, psi, dcosa, toa, spotshape';
for i = 1:numbands
    cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
end
cmd = [cmd,');'];
disp(cmd)
eval(cmd);
                   

disp(Fspot)

