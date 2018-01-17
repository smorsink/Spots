disp('Output results to file!')

cd '/home/kitung/Spot/Data/Data-2017-11-24-X'
%cd '/home/kitung/Spot/FerretData'

load('OptimalSolutions.mat');
%load('MergedSolutions.mat');


para = transpose(OptimalSolutions.X);
dlmwrite('optimals.txt',para,'delimiter', '\t');

fitness = transpose(OptimalSolutions.F);
dlmwrite('fitness.txt',fitness,'delimiter', '\t');

rank = transpose(OptimalSolutions.rank);
dlmwrite('rank.txt',rank,'precision',6);

ndim = size(para)

dlmwrite('dimensions.txt',ndim,'delimiter','\t','precision',6);


