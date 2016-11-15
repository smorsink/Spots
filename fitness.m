function [Fspot, auxOutput] = fitness(X,extPar)
auxOutput = cell(1,size(X,2));
%auxOutput = cell(1,size(X,5));
for i=size(X,2):-1:1  % we want to count backwards here
    % Apply constraints
    isPhysical=applyConstraints(X(:,i),extPar);
    disp('------ Fitness ------');
    if isPhysical
        % Parameters from vector X.
        radius=X(1,i);
        mass=X(2,i);
        inclination=X(3,i);
        emission=X(4,i);
        rho=extPar.fixed.rho;
        timeShift=X(5,i);
        %radius=extPar.fixed.radius;
        %mass=extPar.fixed.mass;
        %inclination=extPar.fixed.inclination;
        %emission=extPar.fixed.emission;
        %timeShift=X(1,i);
        
% Call Mex function in order of fixed quantities
        [Fspot(i),auxOutput{i}] = spotMex_trial(mass, radius, extPar.fixed.freq, inclination, emission, timeShift, ...
            extPar.obsdata2.numbins, extPar.modelchoice, rho, extPar.fixed.spot_temperature, extPar.fixed.distance, ...
            extPar.fixed.numtheta, extPar.spectral_model, extPar.numbands, extPar.fixed.E_band_lower_1, ...
            extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.spots_2, extPar.obsdata2.t, extPar.obsdata2.f(1,:), ...
            extPar.obsdata2.f(2,:), extPar.obsdata2.err(1,:), extPar.obsdata2.err(2,:)); 
    else
        Fspot(i)=Inf;
        auxOutput{i}=[];
        % disp('not physical');
        % disp('');
    end
    
    pause(0.001);
end

function isPhysical=applyConstraints(X, extPar)

parameterNames=extPar.QubistPar.general.XLabels;
m=X(strcmpi(parameterNames,'mass (M_{sun})'));  % the label here needs to match the one in FerretSetup
r=X(strcmpi(parameterNames,'radius (km)'));     % the label here needs to match the one in FerretSetup
%mr_ratio=extPar.fixed.mass/extPar.fixed.radius;
mr_ratio=m/r;
isPhysical=mr_ratio > 0.0677 && mr_ratio < 0.203;