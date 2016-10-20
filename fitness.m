function [Fspot, auxOutput] = fitness(X,extPar)
auxOutput = cell(1,size(X,2));
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
        %inclination=X(1,i);
        %emission=X(2,i);
        %timeShift=X(3,i);
        
% Call Mex function in order of fixed quantities
        [Fspot(i),auxOutput{i}] = spotMex(extPar.obsdata2.numbins, extPar.obsdata2.t, ...
            extPar.modelchoice, mass, radius, extPar.fixed.freq, inclination, ...
            emission, rho, extPar.fixed.spot_temperature, extPar.fixed.E_band_lower_1, ...
            extPar.fixed.E_band_upper_1, extPar.fixed.E_band_lower_2, ...
            extPar.fixed.E_band_upper_2, extPar.fixed.anisotropy, timeShift, ...
            extPar.fixed.bbrat, extPar.fixed.gray, extPar.fixed.Gamma1, ...
            extPar.fixed.Gamma2, extPar.fixed.Gamma3, extPar.fixed.distance, ...
            extPar.obsdata2.f(1,:), extPar.obsdata2.f(2,:), extPar.obsdata2.err(1,:), ...
            extPar.obsdata2.err(2,:));
    else
        Fspot(i)=Inf;
        auxOutput{i}=[];
        disp('not physical');
        disp('');
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