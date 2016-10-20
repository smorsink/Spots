function postProcessing(OptimalSolutions)
disp('------ postProcessing ------');
%global myResults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THIS IS WHAT IS AT THE FRONT OF THE FILE NAMES IN ALL SUBSEQUENT THINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Naming convention: Run + outdata file number + 'a' for first run with that outdata, 'b' for second, etc.
filenameheader='Runx';
% Need to also change this in init.m, outputFerret.m, and spotMex.cpp

exe_dir='/Users/jasper/Documents/spot';

run_shell=strcat(exe_dir,'/contours/',filenameheader,'_ferret_done.sh');
disp(run_shell);
unix(run_shell); % executes the shell script ferret_done.sh


% the rest here writes and executes a short script to clean up the FerretData folder,
%  by cataloguing the files made from the most recent run into a folder
%   with the name 'filenameheader'
temp_str=strcat(exe_dir,'/FerretData/',filenameheader);

cmd1=['mkdir ',temp_str]; % need to do this in order to keep the space after mkdir
disp(cmd1);
unix(cmd1);

%str0=strcat('cd ',exe_dir,'/FerretData \n');
%str1=strcat('mv extPar.mat ',filenameheader,'/extPar.mat \n');
%str2=strcat('mv FerretSetup_SAVED.m ',filenameheader,'/FerretSetup_SAVED.m \n');
%str3=strcat('mv History ',filenameheader,'/History \n');
%str4=strcat('mv init_SAVED.m ',filenameheader,'/init_SAVED.m \n');
%str5=strcat('mv movieFrames ',filenameheader,'/movieFrames \n');
%str6=strcat('mv Notes.txt ',filenameheader,'/Notes.txt \n');
%str7=strcat('mv OptimalSolutions.mat ',filenameheader,'/OptimalSolutions.mat \n');
%str8=strcat('mv par.mat ',filenameheader,'/par.mat \n');
%str9=strcat('mv saveData ',filenameheader,'/saveData \n');
%str10=strcat('mv scratch ',filenameheader,'/scratch \n');
%str11=strcat('mv startTime.mat ',filenameheader,'/startTime.mat \n');
%str12=strcat('echo "FerretData cleaned up. Qubist run files moved to ',filenameheader,'" \n');

movefiles=strcat(exe_dir,'/FerretData/cleanupdirectory.sh');
fileID1=fopen(movefiles, 'w+');
fprintf(fileID1, '# Written in postProcessing.m\n');
fprintf(fileID1, '# Moves all the files and folders created in the most recent Qubist run\n');
fprintf(fileID1, '#  into the appropriate filenameheader folder.\n');
fprintf(fileID1, '# \n');
fprintf(fileID1, 'cd %s/FerretData \n', exe_dir);
fprintf(fileID1, 'mv extPar.mat %s/extPar.mat \n', filenameheader);
fprintf(fileID1, 'mv FerretSetup_SAVED.m %s/FerretSetup_SAVED.m \n', filenameheader);
fprintf(fileID1, 'mv History %s/History \n', filenameheader);
fprintf(fileID1, 'mv init_SAVED.m %s/init_SAVED.m \n', filenameheader);
fprintf(fileID1, 'mv movieFrames %s/movieFrames \n', filenameheader);
fprintf(fileID1, 'mv Notes.txt %s/Notes.txt \n', filenameheader);
fprintf(fileID1, 'mv OptimalSolutions.mat %s/OptimalSolutions.mat \n', filenameheader);
fprintf(fileID1, 'mv par.mat %s/par.mat \n', filenameheader);
fprintf(fileID1, 'mv saveData %s/saveData \n', filenameheader);
fprintf(fileID1, 'mv scratch %s/scratch \n', filenameheader);
fprintf(fileID1, 'mv startTime.mat %s/startTime.mat \n', filenameheader);
fprintf(fileID1, 'echo "FerretData cleaned up. Qubist run files moved to FerretData/%s folder." \n', filenameheader);
fclose(fileID1);

cmd2=strcat('chmod u+rwx "',movefiles,'"');
unix(cmd2);
unix(movefiles);

cmd3=strcat('cp "',exe_dir,'/run_data/',filenameheader,'_texTable.txt" "/Users/jasper/Dropbox/Ferret_runs/',filenameheader,'_texTable.txt"');
unix(cmd3);
cmd4=strcat('cp "',exe_dir,'/run_data/',filenameheader,'_outputFerret_best.txt" "/Users/jasper/Dropbox/Ferret_runs/',filenameheader,'_outputFerret_best.txt"');
unix(cmd4);

disp('postProcessing.m is done.');

