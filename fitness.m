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
        timeShift=X(5,i);
        
        %{
        background1=X(6,i);
        background2=X(7,i);
        background3=X(8,i);
        background4=X(9,i);
        background5=X(10,i);
        %}
        
        
        %radius=extPar.fixed.radius;
        %mass=extPar.fixed.mass;
        %inclination=extPar.fixed.inclination;
        %emission=extPar.fixed.emission;
        %timeShift=extPar.fixed.phaseshift;
        
        
        background1=extPar.fixed.background(1);
        background2=extPar.fixed.background(2);
        background3=extPar.fixed.background(3);
        background4=extPar.fixed.background(4);
        background5=extPar.fixed.background(5);
        
% Call Mex function in order of fixed quantities
        %{
        cmd = '[Fspot(i),auxOutput{i}] = spotMex_trial(mass, radius, extPar.fixed.freq, inclination, emission, timeShift, extPar.obsdata2.numbins, extPar.modelchoice, rho, extPar.fixed.spot_temperature, extPar.fixed.distance, extPar.fixed.numtheta, extPar.spectral_model, extPar.fixed.numbands, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.spots_2, extPar.obsdata2.t';
        for j = 1:extPar.fixed.numbands
            cmd = [cmd,', extPar.obsdata2.f(',num2str(j),',:), extPar.obsdata2.err(',num2str(j),',:)'];
        end
        cmd = [cmd,');'];
        disp(cmd);
        eval(cmd);
        %}
        
        cmd = '[Fspot(i),auxOutput{i}] = spotMex_trial(mass, radius, extPar.fixed.freq, inclination, emission, timeShift, extPar.fixed.numbins, extPar.fixed.modelchoice, extPar.fixed.rho, extPar.fixed.spot_temperature, extPar.fixed.distance, extPar.fixed.numtheta, extPar.fixed.spectral_model, extPar.fixed.numbands, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.fixed.spots_2, extPar.obsdata2.t, extPar.fixed.bend_file_is, extPar.fixed.mr, extPar.fixed.b, extPar.fixed.psi, extPar.fixed.dcosa, extPar.fixed.toa, extPar.fixed.spotshape';
        for j = 1:extPar.fixed.numbands
            cmd = [cmd,', extPar.obsdata2.f(',num2str(j),',:), extPar.obsdata2.err(',num2str(j),',:), background',num2str(j),''];
        end
        cmd = [cmd,');'];
        disp(cmd)
        disp(background1);
        eval(cmd);
    
    
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
%isPhysical=mr_ratio > 0.0677 && mr_ratio < 0.203;
isPhysical=  mr_ratio < 0.175;