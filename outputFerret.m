function outputFerret(world)
disp('------ outputFerret ------');
% Runs once at the end of each generation.

%%%% Finding the minimum chisquared (i.e. best fit) for the generation
FMin=Inf; % double, minimum chisquared of the generation
p_min=0;  % int, index of the population that held the min chisquared of the whole generation
i_min=0;  % int, index of the individual that was had the min chisquared of the whole generation, lives in population p_min

for p=1:length(world.pop)
	%FMin=min(FMin, min(world.pop{p}.indiv.F));
	if min(world.pop{p}.indiv.F) < FMin
		p_min=p;
		[FMin,i_min]=min(world.pop{p_min}.indiv.F);
	end
end
XMin=world.pop{p_min}.indiv.genome.XPhys(:,i_min);

if XMin(5) < 0
	ts = 1 + XMin(5);
else
	ts = XMin(5);
end
%radius = world.par.user.extPar.fixed.radius;
%mass = world.par.user.extPar.fixed.mass;
%incl = XMin(1);
%emis = XMin(2);
radius = XMin(1);
mass = XMin(2);
incl = XMin(3);
emis = XMin(4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS IS WHAT IS AT THE FRONT OF THE FILE NAMES IN ALL SUBSEQUENT THINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Naming convention: Run + outdata file number + 'a' for first run with that outdata, 'b' for second, etc.
filenameheader='Runx';
% Need to also change this in init.m, postProcessing.m, and spotMex.cpp
%
%                         NEED TO CHANGE INPUT FILE HERE TOO 
infile='/Users/jasper/Documents/spot/input/Outdata1.txt';
%
%


%%%% Printing the best fit parameters to a file
filename1=strcat('/Users/jasper/Documents/spot/run_data/',filenameheader,'_outputFerret_best.txt');
fileID=fopen(filename1, 'a+');
fprintf(fileID, '%f \t %f \t %f \t %f \t %f \t %f\n', mass, radius, incl, emis, ts, FMin);
fclose(fileID);

data_dir_name=strcat('/Users/jasper/Documents/spot/FerretData/',filenameheader,'/Notes.txt');
%%%% Printing the best fit parameters to the Notes file, to show up in the Ferret window
notes_file=fopen('/Users/jasper/Documents/spot/FerretData/Notes.txt','a+');
%fprintf(notes_file, 'M= %6f \nR = %6f \ni = %6f \ne = %6f \nts = %6f \n', ...
%mass, radius, incl, emis, ts);
fprintf(notes_file, 'i = %6f \ne = %6f \nts = %6f \n', ...
incl, emis, ts);
fprintf(notes_file, 'chisquared = %8f \n\n', FMin);
fclose(notes_file);

%%%% Printing and saving the light curve plots for the best fit parameters

auxMin=world.pop{p_min}.indiv.auxOutput{i_min};
%disp(auxMin);
numbins = world.par.user.extPar.obsdata2.numbins;

data_time = world.par.user.extPar.obsdata2.t;
data_flux_1 = world.par.user.extPar.obsdata2.f(1,:); % low energy band
data_err_1 = world.par.user.extPar.obsdata2.err(1,:); % error in flux
data_flux_2 = world.par.user.extPar.obsdata2.f(2,:); % high energy band
data_err_2 = world.par.user.extPar.obsdata2.err(2,:);

nothing = zeros(1,numbins); % no error in time

dt = transpose(data_time); % makes it a vertical array
df1 = transpose(data_flux_1); % makes it a vertical array
df2 = transpose(data_flux_2); % makes it a vertical array
de1 = transpose(data_err_1); % makes it a vertical array
de2 = transpose(data_err_2); % makes it a vertical array

dataM = horzcat(dt, df1, df2); % making a matrix of the data values
%disp(dataM);

% Getting the best fit from auxMin (printed to in SpotMex.cpp)
fit_time = auxMin(:,1);
fit_flux_1 = auxMin(:,2);
fit_flux_2 = auxMin(:,3);

fitM = horzcat(fit_time, fit_flux_1, fit_flux_2); % making a matrix of the fit values
%disp(fitM);


%% For file naming writing purposes
intmass=cast(mass*1e6, 'uint64');
intradius=cast(radius*1e6, 'uint64');
intincl=cast(incl*1e6, 'uint64');
intemis=cast(emis*1e6, 'uint64');
intts=cast(ts*1e6, 'uint64');


%disp(intmass);
%disp(intradius);
%disp(intincl);
%disp(intemis);
%disp(intts);
spot_output_file=strcat('/Users/jasper/Documents/spot/output/ferret/q',int2str(world.par.user.extPar.modelchoice),'_g',int2str(world.par.user.extPar.fixed.gray),'_m',int2str(intmass),'e-6_r',int2str(intradius),'e-6_i',int2str(intincl),'e-6_e',int2str(intemis),'e-6_ts',int2str(intts),'e-6.txt');
%disp(spot_output_file);


%filename3=strcat('run_data/',filenameheader,'_outputFerret_fit.txt');
%filename4=strcat('run_data/',filenameheader,'_outputFerret_data.txt');
filename5=strcat('run_data/',filenameheader,'_ferret_to_spot_cmdline.txt');
% Writing matrices to output files, for checking values
% not working when I want to put a few header lines in
%fitprint=fopen(filename3,'w+');
%fprintf(fitprint,'# Fit values from myPlot \n');
%fprintf(fitprint,'# Column 1: time \n# Column 2: Flux, first energy band \n# Column 3: Flux, second energy band \n');
%fprintf(fitprint,'# \n');
%fclose(fitprint);
%dlmwrite(filename3, fitM, '\t');
%dataprint=fopen(filename4,'w+');
%fprintf(dataprint,'# Data read into myPlot \n');
%fprintf(dataprint,'# Column 1: time \n# Column 2: Flux, first energy band \n# Column 3: Flux, second energy band \n');
%fprintf(dataprint,'# \n');
%fclose(dataprint);
%dlmwrite(filename4, dataM, '\t');

%% writing ferret_to_spot_cmdline.txt
fileID5=fopen(filename5, 'a+');
fprintf(fileID5,'"/Users/jasper/Documents/spot/spot" -m %f -r %f -f %d -i %f -e %f -l %f -n 32 -q %u -o "%s" -p %g -T %g -D %g -t 1 -u %g -U %g -v %g -V %g -g %u -I "%s" -3 "%s" -N \n', ...
mass, radius, world.par.user.extPar.fixed.freq, incl, emis, ts, world.par.user.extPar.modelchoice, spot_output_file, world.par.user.extPar.fixed.rho, world.par.user.extPar.fixed.spot_temperature, world.par.user.extPar.fixed.distance, ...
world.par.user.extPar.fixed.E_band_lower_1, world.par.user.extPar.fixed.E_band_upper_1, world.par.user.extPar.fixed.E_band_lower_2, world.par.user.extPar.fixed.E_band_upper_2, world.par.user.extPar.fixed.gray, infile, filenameheader);
fclose(fileID5);

%% FOR PARAM_DEGEN
PD_div=50;
m_min=1.0;
m_max=2.5;
r_min=10.0;
r_max=16.0;

%% writing Ferret_to_pd_cmdline.txt
filename6=strcat('/Users/jasper/Documents/spot/contours/',filenameheader,'_ferret_to_pd_cmdline.txt');
fileID6=fopen(filename6, 'a+');
fprintf(fileID6,'"/Users/jasper/Documents/spot/param_degen/param_degen" -a %g -A %g -b %g -B %g -c %f -m %f -r %f -f %d -i %f -e %f -l %f -n 32 -q %u -o "%s" -p %g -T %g -D %g -u %g -U %g -v %g -V %g -g %u -I "%s" -3 "%s" -d %u \n', ...
m_min, m_max, r_min, r_max, FMin, mass, radius, world.par.user.extPar.fixed.freq, incl, emis, ts, world.par.user.extPar.modelchoice, spot_output_file, world.par.user.extPar.fixed.rho, world.par.user.extPar.fixed.spot_temperature, world.par.user.extPar.fixed.distance, ...
world.par.user.extPar.fixed.E_band_lower_1, world.par.user.extPar.fixed.E_band_upper_1, world.par.user.extPar.fixed.E_band_lower_2, world.par.user.extPar.fixed.E_band_upper_2, world.par.user.extPar.fixed.gray, infile, filenameheader, PD_div);
fclose(fileID6);

low_from_mex=strcat('/Users/jasper/Documents/spot/contours/',filenameheader,'_low_from_spotmex.txt'); % don't forget to change this in spotMex too
spot_for_gnuplot=strcat('/Users/jasper/Documents/spot/contours/',filenameheader,'_spot_for_gnuplot.txt'); % don't forget to change this in paramdegen and spot


%% writing Ferret_done.sh
filename7=strcat('/Users/jasper/Documents/spot/contours/',filenameheader,'_ferret_done.sh');
fileID7=fopen(filename7, 'w+');
fprintf(fileID7, '# Written to in outputFerret.m\n');
fprintf(fileID7, '# To be run at the end of a full Ferret run \n');
fprintf(fileID7, '# Overwritten every time outputFerret runs, so that it is left with just the last run of outputFerret.\n');
fprintf(fileID7, '# \n');
fprintf(fileID7, 'cd /Users/jasper/Documents/spot/param_degen\n');
fprintf(fileID7, 'make param_degen\n');
fprintf(fileID7,'"/Users/jasper/Documents/spot/param_degen/param_degen" -a %g -A %g -b %g -B %g -c %f -m %f -r %f -f %d -i %f -e %f -l %f -n 32 -q %u -o "%s" -p %g -T %g -D %g -u %g -U %g -v %g -V %g -g %u -I "%s" -3 "%s" -d %u \n', ...
m_min, m_max, r_min, r_max, FMin, mass, radius, world.par.user.extPar.fixed.freq, incl, emis, ts, world.par.user.extPar.modelchoice, spot_output_file, world.par.user.extPar.fixed.rho, world.par.user.extPar.fixed.spot_temperature, world.par.user.extPar.fixed.distance, ...
world.par.user.extPar.fixed.E_band_lower_1, world.par.user.extPar.fixed.E_band_upper_1, world.par.user.extPar.fixed.E_band_lower_2, world.par.user.extPar.fixed.E_band_upper_2, world.par.user.extPar.fixed.gray, infile, filenameheader, PD_div);
fprintf(fileID7, 'cd /Users/jasper/Documents/spot/chi2contours\n');
fprintf(fileID7, 'make chi2contours\n');
fprintf(fileID7, '"/Users/jasper/Documents/spot/chi2contours/chi2contours" -c %f -r %f -m %f -3 "%s" -i "%s" -I "%s" \n', FMin, radius, mass, filenameheader, low_from_mex, spot_for_gnuplot);
fprintf(fileID7, 'cd /Users/jasper/Documents/spot/gridded\n');
fprintf(fileID7, 'make griddedchi2\n');
fprintf(fileID7, '"/Users/jasper/Documents/spot/gridded/griddedchi2" -a %g -A %g -b %g -B %g -3 "%s" -i "%s" -d %u \n',m_min, m_max, r_min, r_max, filenameheader, low_from_mex, PD_div);
fclose(fileID7);

filename8=strcat('/Users/jasper/Documents/spot/run_data/',filenameheader,'_texTable.txt');
fileID8=fopen(filename8, 'w+');
fprintf(fileID8,' & Fit %s & %.2f & %.2f & %.2f & %.2f & %.2f & %.4f & \\\\  %% %s \n', filenameheader, FMin, mass, radius, incl, emis, ts, filenameheader);
fprintf(fileID8,' & Fit %s & %.2f & %.2f & \\%% & %.2f & \\%% & %.2f & %.2f & %.4f \\\\  %% %s \n', filenameheader, FMin, mass, radius, incl, emis, ts, filenameheader);
fclose(fileID8);

% Setting up the figure for saving
filename2=strcat('run_data/',filenameheader,'_fitplot.ps');
h=figure;

errorbar(dt,df2,nothing,de2);
hold on
plot(fit_time,fit_flux_2,'-m');

print(h,'-dpsc2','-append',filename2); % this has a plot for every generation of the run

%saveas(h,'-psc2',filename2); % this one only has the most recent one

% Unique stopping criteria
if FMin <= 59.0
	abortQubist(1);
end
