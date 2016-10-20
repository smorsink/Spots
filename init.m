function extPar = init  %defining an external parameter sent to the executable; equivalent of command line arguments to main

spotDir=fileparts(which(mfilename));

% ------------
% Compile the functions...
%
% If you see bizarre errors when mex evaluates..
% There is a file called mexopts.sh that needs to be modified.
% It''s located in /home/<your user name>/.matlab/R2010a/mexopts.sh (may depend on MATLAB version)
% look in there and you'll see about four places where "-ansi" is specified as a flag.
% Simply remove the "-ansi" and it should compile.
%
cmd=['mex spotMex.cpp -L', spotDir, ' -cxx Units.cpp -cxx Chi.cpp -cxx OblDeflectionTOA.cpp -cxx PolyOblModelBase.cpp -cxx PolyOblModelCFLQS.cpp -cxx PolyOblModelNHQS.cpp -cxx matpack.cpp -cxx SphericalOblModel.cpp'];

disp('------ Init ------');
disp(cmd)
which('spotMex.cpp');

eval(cmd);
%
%

% ------------
% extPar is external parameters
%
extPar.modelchoice=3
;               % 1 is for oblate, 3 is for spherical
%
extPar.fixed.freq=581;              % spin frequency of NS, in Hz
%
extPar.fixed.E_band_lower_1=2.0;    % lower bound of first energy band, in keV
extPar.fixed.E_band_upper_1=3.0;    % upper bound of first energy band, in keV
extPar.fixed.E_band_lower_2=5.0;    % lower bound of second energy band, in keV
extPar.fixed.E_band_upper_2=6.0;    % upper bound of second energy band, in keV
extPar.fixed.bbrat=1;               % ratio of blackbody to comptonization (default = 1)
extPar.fixed.anisotropy=0.586;      % anisotropy parameter (default = 0.586)
extPar.fixed.spot_temperature=2.0;  % temperature of the hot spot in the NS's frame, in keV
extPar.fixed.distance=1.85141655e20;% distance from us to the NS, in meters
extPar.fixed.gray=0;                % toggle graybody factor in spectrum, 0=no, 1=yes (default = 0)
extPar.fixed.rho=6.0;               % angular radius of the emitting hot spot, in degrees
%
%extPar.fixed.mass=1.465;
%extPar.fixed.radius=11.68;

%
extPar.fixed.Gamma1=2.0;
extPar.fixed.Gamma2=2.0;
extPar.fixed.Gamma3=2.0;
%
% ------------
% Load data.
%               NEED TO CHANGE THIS IN outputFerret.m TOO, ~ line 80
exe_dir='/Users/jasper/Documents/spot';
fileName=strcat(exe_dir,'/input/Outdata1.txt'); % *** need to change this in outputFerret too
dataFile=load(fileName);
% ------------
% Process data
% For file 1998_box4
%extPar.obsdata1.t=dataFile(:,1)';
%extPar.obsdata1.f(1,:)=dataFile(:,4)';
%extPar.obsdata1.err(1,:)=dataFile(:,5)';
%extPar.obsdata1.f(2,:)=dataFile(:,8)';
%extPar.obsdata1.err(2,:)=dataFile(:,9)';
%extPar.obsdata1.numbins=size(dataFile,1);
%
% For file test_input_#
extPar.obsdata2.t=dataFile(:,1)';
extPar.obsdata2.f(1,:)=dataFile(:,2)';
extPar.obsdata2.err(1,:)=dataFile(:,3)';
extPar.obsdata2.f(2,:)=dataFile(:,4)';
extPar.obsdata2.err(2,:)=dataFile(:,5)';
extPar.obsdata2.numbins=size(dataFile,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS IS WHAT IS AT THE FRONT OF THE FILE NAMES IN ALL SUBSEQUENT THINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Naming convention: Run + outdata file number + 'a' for first run with that outdata, 'b' for second, etc.
filenameheader='Runx';
%% need to also change this in outputFerret.m, postProcessing.m, and spotMex.cpp


%% need to have full filenames here, otherwise get fopen error in matlab

filename2=strcat(exe_dir,'/run_data/',filenameheader,'_ferret_to_spot_cmdline.txt');
fileID2=fopen(filename2, 'a+');
fprintf(fileID2, '# Written in init.m and outputFerret.m\n');
fprintf(fileID2, '# Input file: %s \n', fileName);
fprintf(fileID2, '# \n');
fclose(fileID2);

% writes header for everything_from_spotmex.txt
filename3=strcat(exe_dir,'/contours/',filenameheader,'_everything_from_spotmex.txt'); % don't forget to change this in spotMex too
fileID3=fopen(filename3, 'a+');
fprintf(fileID3,'# Constants: f = %d Hz, T = %g keV, p = %g degrees, D = %g m, gray = %u\n', ...
extPar.fixed.freq, extPar.fixed.spot_temperature, extPar.fixed.rho, extPar.fixed.distance, extPar.fixed.gray);
fprintf(fileID3, '# \t\tNS model = %u, low energy band = %g-%g keV, high energy band = %g-%g keV\n', ...
extPar.modelchoice, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.E_band_lower_2, extPar.fixed.E_band_upper_2);
fprintf(fileID3, '# \t\tInput file = %s\n', fileName);
fprintf(fileID3, '# Column 1: Mass (Msun) \n');
fprintf(fileID3, '# Column 2: Radius (km) \n');
fprintf(fileID3, '# Column 3: Inclination angle (degrees) \n');
fprintf(fileID3, '# Column 4: Emission angle (degrees) \n');
fprintf(fileID3, '# Column 5: Phase shift \n');  
fprintf(fileID3, '# Column 6: Chisquared \n');
fprintf(fileID3, '# \n');
fclose(fileID3);

% writes header for lowchisquared_from_spotmex.txt
filename4=strcat(exe_dir,'/contours/',filenameheader,'_low_from_spotmex.txt'); % don't forget to change this in spotMex and outputferret
fileID4=fopen(filename4, 'a+');
fprintf(fileID4,'# Constants: f = %d Hz, T = %g keV, p = %g degrees, D = %g m, gray = %u\n', ...
extPar.fixed.freq, extPar.fixed.spot_temperature, extPar.fixed.rho, extPar.fixed.distance, extPar.fixed.gray);
fprintf(fileID4, '# \t\tNS model = %u, low energy band = %g-%g keV, high energy band = %g-%g keV\n', ...
extPar.modelchoice, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.E_band_lower_2, extPar.fixed.E_band_upper_2);
fprintf(fileID4, '# \t\tInput file = %s\n', fileName);
fprintf(fileID4, '# Column 1: Mass (Msun) \n');
fprintf(fileID4, '# Column 2: Radius (km) \n');
fprintf(fileID4, '# Column 3: Inclination angle (degrees) \n');
fprintf(fileID4, '# Column 4: Emission angle (degrees) \n');
fprintf(fileID4, '# Column 5: Phase shift \n');  
fprintf(fileID4, '# Column 6: Chisquared \n');
fprintf(fileID4, '# \n');
fclose(fileID4);

% write header for ferret_to_paramdegen.sh
filename5=strcat(exe_dir,'/contours/',filenameheader,'_ferret_to_pd_cmdline.txt');
fileID5=fopen(filename5, 'a+');
fprintf(fileID5, '# Written in init.m and outputFerret.m\n');
fprintf(fileID5, '# \n');
fclose(fileID5);

% write header for paramdegen_to_spot.sh
filename6=strcat(exe_dir,'/contours/',filenameheader,'_pd_to_spot.sh'); % don't forget to change this in paramdegen too
fileID6=fopen(filename6, 'a+');
chmodit6=strcat('chmod u+rwx "',filename6,'"'); % writes the unix command to chmod the file for unix
unix(chmodit6); % executes the chmod command
fprintf(fileID6, '# Written in init.m and param_degen.cpp\n');
fprintf(fileID6, '# Ferret input file: %s \n', fileName);
fclose(fileID6);

% writes header for outputFerret_best.txt
filename7=strcat(exe_dir,'/run_data/',filenameheader,'_outputFerret_best.txt');
fileID7=fopen(filename7,'a+');
fprintf(fileID7, '# Printed to from init.m and outputFerret.m \n');
fprintf(fileID7, '# Parameters with best chisquared of each generation. \n');
fprintf(fileID7, '# Input file = %s\n', fileName);
fprintf(fileID7,'# Constants: f = %d Hz, T = %g keV, p = %g degrees, D = %g m, gray = %u\n', ...
extPar.fixed.freq, extPar.fixed.spot_temperature, extPar.fixed.rho, extPar.fixed.distance, extPar.fixed.gray);
fprintf(fileID7, '# \t\tNS model = %u, low energy band = %g-%g keV, high energy band = %g-%g keV\n', ...
extPar.modelchoice, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.E_band_lower_2, extPar.fixed.E_band_upper_2);
fprintf(fileID7, '# Column 1: Mass (Msun) \n');
fprintf(fileID7, '# Column 2: Radius (km) \n');
fprintf(fileID7, '# Column 3: Inclination angle (degrees) \n');
fprintf(fileID7, '# Column 4: Emission angle (degrees) \n');
fprintf(fileID7, '# Column 5: Phase shift \n');  
fprintf(fileID7, '# Column 6: Chisquared \n');
fprintf(fileID7, '# \n');
fclose(fileID7);

% write header for ferret_done.sh
filename8=strcat(exe_dir,'/contours/',filenameheader,'_ferret_done.sh');
fileID8=fopen(filename8, 'a+');
chmodit8=strcat('chmod u+rwx "',filename8,'"'); % writes the unix command to chmod the file for unix
unix(chmodit8); % executes the chmod command
fprintf(fileID8, '# Written in init.m and outputFerret.m\n');
fprintf(fileID8, '# Ferret input file: %s \n', fileName);
fprintf(fileID8, '# \n');
fclose(fileID8);
