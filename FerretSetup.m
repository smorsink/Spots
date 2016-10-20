% ======================================================================
% Qubist: A Global Optimization, Modeling & Visualization Platform
%
% Ferret: A Multi-Objective Linkage-Learning Genetic Algorithm
% Locust: A Multi-Objective Particle Swarm Optimizer
% Anvil: A Multi-Objective Simulated Annealing/Genetic Algorithm Hybrid
% SAMOSA: Simple Approach to a Multi-Objective Simplex Algorithm
% SemiGloSS: An Experimental Semi-Global Optimizer & Polisher
%
% Copyright 2002-2010 Jason D. Fiege, Ph.D.  All rights reserved.
% Distributed by nQube Technical Computing Corporation under license.
% design.innovate.optimize @ www.nQube.ca
% ======================================================================

function par=FerretSetup(par)
%
% -----------------------------------
% Sets reasonable defaults for parameters where possible.
% Values are over-ridden by FerretSetup.m
% 
% ====================================
% User
disp('------ FerretSetup ------');
par.user.fitnessFcn='fitness'; % [string]: Name of the fitness function.
par.user.output='outputFerret'; % [string]: Name of optional user-defined output function called each generation.
par.user.XMod=''; % [string]: Name of optional user-defined XMod function.
%
% ====================================
% History
%data_dir_name=strcat('FerretData/',date); % the desired name of the history directory
%disp(data_dir_name); % displaying the history directory name to double check that it says what I want it to
% need to change data_dir_name in outputFerret too
par.history.dataDir='FerretData'; % [string]: Directory for Ferret History files, etc.
par.history.NGenPerHistoryFile=25; % [integer >= 1]: How many generations per History file?
par.history.saveFrac=0; % [0 - 1]: Save only the optimals --> [0], or a fraction [0 - 1].  See also par.analysis.keepFrac.
%
% ====================================
% General
% multiple smaller populations is better than one very large population
par.general.NPop=2; % p [integer >= 1]: Number of populations for one generation.
par.general.popSize=3; % i [integer >= 1]: Size of each population.
% for debugging purposes it helps to have smaller numbers for popSize and NPop=1
par.general.NGen=4; % g [integer >= 1]: Maximum number of generations to run for.
% pop=15, indiv=300, gen=30
% OR pop=20, indiv=350, gen=35
par.general.FLabels={'\chi^2'}; % [Cell array of strings]: Give names to some or all fitness values: {'FA','FB',...}
%
% -m mass in Msun
% -r radius (at spot) in km
% -f frequency in Hz
% -i inclination angle (degrees)
% -e emission angle, angle from geographic north to center of hot spot (degrees)
% -p angular radius of spot (degrees)
% -T temperature of the spot, in the spot''s frame (keV)
% -a anisotropy parameter
% -l phase shift (between 0 and 1)
% -b bbrat
%
par.general.XLabels={'radius (km)', 'mass (M_{sun})', 'inclination (degrees)', '\theta (degrees)', 'phase shift'};
par.general.min=[       10.0,           1.0,                 0.01,                 0.01,               0.00     ];
par.general.max=[       16.0,           2.5,                 90.0,                 90.0,               1.00     ];
%r:10-16, m:1-2.5, i:0.01-90, e:0.01-90, ts:0-1
%par.general.XLabels={'inclination (degrees)', '\theta (degrees)', 'phase shift'};
%par.general.min=[          0.01,                  0.01,               0.00     ];
%par.general.max=[          90.0,                  90.0,               1.00     ];
par.general.softMin=[]; % [integer vector > 1]: Which parameters have soft minimum bounds?
par.general.softMax=[]; % [integer vector > 1]: Which parameters have soft maximum bounds?
par.general.discreteStep=[]; % [real vector >= 0]: Discrete stepsizes of some parameters?  0's or [] --> continuous.
par.general.cyclic=[5]; % [integer vector > 1]: Which parameters are cyclic?
%par.general.cyclic=[3]; % [integer vector > 1]: Which parameters are cyclic?
par.general.HolderMetricExponent=2; % [real]: Holder metric for determining distances.
par.general.noisy=false; % [logical]: Does the fitness function contain noisy fluctuations?
%
% ====================================
% Strategy Parameters
par.strategy.isAdaptive=true; % [logical]: Global switch for strategy adaptation.
par.strategy.NGen=25; % [integer > 1]: Number of generations over which trackers are smoothed and forget their history.
%
% Which strategy parameters are allowed to adapt? (true or false)
par.strategy.adapt.mutation.scale=true; % [logical]
par.strategy.adapt.XOver.dispersion=true; % [logical]
par.strategy.adapt.XOver.strength=true; % [logical]
par.strategy.adapt.XOver.ALS=true; % [logical]
par.strategy.adapt.XOver.matingRestriction=true; % [logical]
par.strategy.adapt.niching.acceleration=true; % [logical]
%
% ====================================
% Parallel Computing
par.parallel.NWorkers=0; % [integer >= 0]: Number of worker nodes to launch initially.
par.parallel.nodeDistributionFactor=5; % [integer >= 1]: average number of chunks each node evaluates.
par.parallel.minChunkSize=10; % [integer]: Minimum number of evaluations in each work chunk.
par.parallel.timeout=30; % [integer >= 0]: Maximum time in seconds before an unresponsive node disconnects.
par.parallel.latency=0; % [real > 0]: Time to pause in seconds after writing to the scratch directory.
par.parallel.useJava=true; % [logical]: Is Java required for worker nodes?
par.parallel.writeLogFiles=true; % [logical]: Are log files required?
%
% ====================================
% Selection
par.selection.PCompete=1; % [0 - 1]: Probability that each individual will compete.
par.selection.pressure=0.8; % [0 - 1]: Selection pressure on overall fitness.
par.selection.BBPressure=1; % [0 - 1]: Selection pressure on BBs.
par.selection.FAbsTol=0; % Absolute fitness range (+/- FAbsTol) to use as a fuzzy fitness band.
par.selection.FRelTol=0.05; % [0 - 1]: Fraction of fitness range to use as a fuzzy fitness band.
par.selection.FRankTol=0; % [0 - 1]: Fraction of rank range to use as a fuzzy fitness band.
par.selection.exploitFrac=0.10; % [0 - 1]: Population fraction devoted to exploiting the optimal region rather than exploring.
par.selection.shareFitness=0; % [0 - 1]: Fitness sharing.
par.selection.clusterScale=1; % [0 - 1]: Setting clusterScale < 1 promotes clusters on the requested scale.
par.selection.competitionRestriction=0; % [-1:1]: Negative values enhance probability of competitions between nearby individuals.
par.selection.minimizeEachObjective=0; % [real vector]: Probablilities to enhance minimization of each objective separately.
% Multi-objective functions only.  0 --> no extra effort.  1 --> max.
%
% ====================================
% Mutation
par.mutation.PMutate=0.15; % [0 - 1]: Probability of normal mutation.
par.mutation.BBRestricted=false; % [logical]: Restrict mutation to building blocks?
par.mutation.PRandomizeBBs=0.01; % [0 - 1]: Probability of randomizing whole building blocks.
par.mutation.PSuperMutate=0.01; % [0 - 1]: Probability of superMutation.
par.mutation.scale=0.4; % [0 - 1]: Min & max scale of mutation.
par.mutation.selectiveMutation=-0.2; % [-1:1]: Negative values mutate highly clustered individuals selectively.
%
% ====================================
% Crossover
par.XOver.PXOver=1;  % [0 - 1]: Probability of X-Type crossovers.
%previously 1
par.XOver.BBRestricted=false; % [logical]: Restrict crossover to building blocks?
par.XOver.strength=0.5; % [0 - 1]: Scale of crossover.
par.XOver.maxScale=1.1; % [~1, or slightly larger]: maximum XOver perturbation allowed.
par.XOver.antiXOverProb=0; % [0 - 1]: Probability of anti-crossovers.
par.XOver.ALS=0; % [0 - 1]: ALS: Advanced Lethal Suppression.
% ------------------------------------------
% Which dispersion technique should be used?
% par.XOver.dispersionTechnique='cylindrical'; % [string]: Hollow cylinder.
% par.XOver.dispersionTechnique='conic'; % [string]: Scaled hollow cone.
par.XOver.dispersionTechnique='biconic'; % [string]: Scaled hollow double-cone.
% ------------------------------------------
par.XOver.dispersion=0.25; % [real < 1]: Determines major axis standard deviation of dispersion structure.
par.XOver.matingRestriction=0; % [-1:1]: Negative values enhance probability of matings between
% relatively nearby individuals.  Positive values are probably not useful.
%
% ====================================
% Building Block Crossover
par.XOverBB.PXOver=1; % [0 - 1]: Probability of Building Block crossovers.
par.XOverBB.multiBB=false; % [logical]: Single or multiple BB XOver?
par.XOverBB.NPass=1; % [integer >= 0]: Number of passes through BB selection per cycle.
%
% ====================================
% Niching
par.niching.priority='PXF'; % [string: re-order 3 letters, 'P', 'X', or 'F']: Priority of niching. If empty, use Pareto niching.
par.niching.P=0.5;  % [0 - 1]: Pattern niching: Typically ~0.5 is about right.
par.niching.X=0.25;  % [0 - 1]: X-Niching: Typically ~0.25 is about right.
par.niching.XPar=[]; % [integer vector > 1]: List of parameters used for X-niching.  If empty, use all.
par.niching.F=0.25;  % [0 - 1]: F-Niching: Typically ~0.25 is about right.
par.niching.method='sigmaShare'; % [string: 'sigmaShare' or 'powerLaw']: Specify a niching method.
par.niching.exponent=2; % [real: usually > 0 & <~ 2]: Used in the niche function.
par.niching.powers=[4]; % [real vector: usually > 0]: Used only for power-law niching - not sigmaShare. 
par.niching.acceleration=0.5; % [0 - 1]: Acceleration parameter.
%
% ------------------------------------
% Pattern Niching Specifics
par.niching.patternMethod='sigmaShare'; % [string: 'sigmaShare' or 'powerLaw']: Method to use for pattern niching.
par.niching.patternExponent=2; % [real: usually > 0 & <~ 2]: exponent for pattern niching.
par.niching.patternPowers=[2]; % [real vector, usually > 0]: Used only for power-law niching.
%
% ====================================
% Critical Parameter Detection
par.CPD.PDeactivateGenes=0.000; % [0 - 1]: Probability of gene deactivation.
par.CPD.PReactivateGenes=0.000; % [0 - 1]: Probability of gene reactivation.
par.CPD.useTemplates=false; % [logical]: Use templates to replace NaN's in evaluation when possible.
par.CPD.NaN2Random=true; %  [logical, usually true]: Replace undefined genes (NaN's) with random numbers?
%
% ====================================
% Immigration
par.immigration.PImmigrate=0.01; % [0 - 1, usually << 1]: Probability of immigration between populations.
%
% ====================================
% Elitism
par.elitism.mode='normal'; % ['normal' or 'boundary']: Elitism mode.
par.elitism.frac=0.05; % [0 - 1]: Fraction of population that are designated elite.
%
% ====================================
% Linkage Learning
par.link.PLink=1; % [0 - 1]: Probability that an individual will be chosen for linkage learning.
par.link.FAbsTol=1; % Absolute fitness range (+/- FAbsTol) to use as a fuzzy fitness band.
par.link.FRelTol=0; % [0 - 1]: Fraction of fitness range to use as a fuzzy fitness band.
par.link.FRankTol=0; % [0 - 1]: Fraction of rank range to use as a fuzzy fitness band.
par.link.maxStrongLinkFrac=1; % [0 - 1]
par.link.NLinkMax=Inf; % [integer >= 1]: Maximum number of linkage attempts per generation.
par.link.lifetime=100; % [integer >= 1]: Lifetime of linkages.
par.link.useOptimalDonors=false; % [logical]: Should a member of the Optimal set be used for the donor?
par.link.acceleration=2; % [real >= 1]: Accelerate decay of weak links.
par.link.failureCost=1; % [real ~1]: Cost associated with higher-order links.
par.link.initialLinkage=0.5; % [0 - 1]: How linked are the parameters initially?
par.link.maxPass=20; % [integer >= 1]: Max number of passes through the linkage queue.
par.link.searchThreshold=0.9; % [0 - 1]: Stop looking for linkages when links are stronger than this value.
par.link.PMultiLink=0.5; % [0 - 1]: Linkage strategy.  0 --> simple links only.  1 --> general, higher order links possible.
par.link.PXOverDonor=1;  % [0 - 1]: Probability to engage in XOver with the donor.
par.link.strategy=0.5; % [0 - 1]: 0 --> opportunistic, 1 --> aggressive.
par.link.PRandomizeNewBBs=0;  % [0 - 1]: Randomize BBs when they form?
par.link.mixDim=0; % [usually 0 - 1, possibly higher]: Greater than 0 if mixing is desired.
par.link.BB={}; % [linkage matrix (0-1), or cell array of integer vectors]: Known building blocks. Empty list if none are known.
%
% ====================================
% Zooming % WE DO NOT WANT ZOOMING!
par.zoom.NGen=0; % [integer >= 1]: Zoom control.
par.zoom.safety=0.5; % [real, usually <= 1]: Size of buffer zone around optimal set.
par.zoom.allowZoomOut=true; % [logical]: Is zoom-out permitted?
par.zoom.maxPower=1e10; % [real << 1]: Maximum zoom power.
%
% ====================================
% Stopping criteria
%
% Tolerance on X and F - mainly for single-objective problems. 
par.stopping.XTol=NaN; % [real > 0, or real vector (length=% of parameters)]: Tolerance on parameters.
par.stopping.FTol=NaN; % [real > 0, or real vector (length=% of objectives)]: Tolerance on objectives.
par.stopping.boolean=@or; % [boolean operator]: Should be '@and' or '@or'.
par.stopping.window=25; % [real ~ 1]: Ferret keeps par.tracks.window generations in memory for stopping criteria.
%
% ====================================
% Analysis
par.analysis.analyzeWhenDone=true; % [logical]: Do analysis automatically when evolution stops.
par.analysis.conserveMemory=false; % [logical]: Minimize memory usage during analysis?  (Analysis will be slower)
par.analysis.maxItNoProgress=250; % [integer >= 1]: Max iterations with no progress before analysis stops.
par.analysis.postProcess='postProcessing'; % [string]: Name of code to process OptimalSolutions when analysis is done.
%
% ====================================
% Local Optimization
par.localOpt.startGen=Inf; % [integer]: First generation for local optimization. (Inf or NaN for no set start generation)
par.localOpt.startF=NaN; % [real]: Trigger local optimization when fitness falls below some value. (NaN for no trigger)
par.localOpt.convergeWindow=50; % [integer > 1]: The length of the window (in generations) used to monitor for convergence.
par.localOpt.startFTol=0; % [real]: Trigger local optimization when the change in fitness drops below this value (single-objective only).
par.localOpt.startWhenConverged=false; % [logical]: Trigger local optimization when population appears converged.
par.localOpt.POptimizeOnImprovement=0; % [0 - 1]: Probability of optimization after improvement of best solution.
par.localOpt.POptimizeOnNoImprovement=0; % [0 - 1]: Probability of optimization after no improvement of best solution.
par.localOpt.maxFEval=1000; % [integer]: Maximum number of evaluations allowed for each local optimization.
par.localOpt.tolX=1e-6; % [real > 0, usually small]: Tolerance in X for local optimization.
par.localOpt.tolF=1e-6; % [real > 0, usually small]: Tolerance in F for local optimization.
% 
% ------------------------------------
% Selection of optimals to undergo local optimization: The most *restrictive* of these 2 criteria are used: i.e the one that yields
% the fewest optimals.  Each time the local optimizer is called, select a fraction par.localOpt.frac of the solutions in the current
% optimal set (i.e. within the fuzzy tolerance of the best solution for single objective problems, or solutions on the Pareto front
% for multi-objective problems) UP TO a maximum number of optimals equal to par.localOpt.maxNOptimals.
%
par.localOpt.frac=NaN; % [0 - 1]: Fraction of optimals to undergo local optimization.  NaN defaults value to a single optimal for
% single-objective problems, and all optimals on the Pareto front for multi-objective problems.  At least 1 optimal is chosen.
par.localOpt.maxNOptimals=NaN; % [integer >= 1]: Maximum number of optimals to undergo local optimization.  NaN or Inf
% means that no additional filtering will be done -- the entire fraction selected (par.localOpt.frac) will be used.
%
% ====================================
% *** Polisher Setup ***
%
% par.polish.optimizer='fminsearch';
% par.polish.optimizer='SAMOSA';
% par.polish.optimizer='Anvil';
% par.polish.optimizer='SemiGloSS';
par.polish.optimizer='default'; % [string]: Use defaults.
%
% fminsearch (Single-objective only):
par.polish.fminsearch.NTracks=Inf; % [integer]: Maximum number of tracks.  `Inf' --> all in optimal set.
par.polish.fminsearch.options=[]; % options structure from optimset.
%
% SAMOSA uses its own setup file.  Here, we just add the number of tracks requested.
par.polish.SAMOSA=defaultSAMOSA_setup; % [m-file name]: Name of SAMOSA setup function.
par.polish.SAMOSA.NTracks=Inf; % [integer]: Maximum number of tracks.  `Inf' --> all in optimal set.
%
% Anvil & SemiGloSS uses their own setup files.
par.polish.Anvil=defaultAnvilSetup;  % [m-file name]: Name of Anvil setup function.
par.polish.SemiGloSS=defaultSemiGloSS_setup; % [m-file name]: Name of SemiGloSS setup function.
%
% ====================================
% Movies
par.movie.makeMovie=false; % [logical]: Generate a movie of the run?
par.movie.frameRate=5; % [integer >= 1]: Movie frame rate per second.
%
% ====================================
% Graphics & messages
par.interface.graphics=1; % [integer: 1, 0, or -1]: Full graphics (1), low graphics (0), or no graphics (-1)?
par.interface.myColorMap='bone'; % [string: colormap name] User choice for colormap.
par.interface.myPlot='myPlot'; % [string]: Name of custom plot function .m file.
par.interface.fontUnits='points'; % [string]: Must be 'points' or 'pixels'.
par.interface.titleFontSize=12; % [integer >= 1]: Self-explanatory.
par.interface.labelFontSize=10; % [integer >= 1]: Self-explanatory.
par.interface.axisFontSize=8; % [integer >= 1]: Self-explanatory.
par.interface.xAxis.type='X'; % ['X' or 'F']: Default X-axis type
par.interface.xAxis.value=1; % [1, NPar] or [1,NObj]: Default X-axis variable
par.interface.yAxis.type='X'; % ['X' or 'F']: Default Y-axis type
par.interface.yAxis.value=2; % [1, NPar] or [1,NObj]: Default Y-axis variable
par.interface.zAxis.type='F'; % ['X' or 'F']: Default Z-axis type
par.interface.zAxis.value=1; % [1, NPar] or [1,NObj]: Default Z-axis variable
par.interface.NPix=25; % [integer >= 1]: Size of the grid for contour and mesh plots.
par.interface.NContours=10; % [integer >= 1]: Number of contour levels for contour plots.
