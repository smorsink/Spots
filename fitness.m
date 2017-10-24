function [Fspot, auxOutput] = fitness(X,extPar)
auxOutput = cell(1,size(X,2));
% This version is for comparison with Slavko's synthetic data
for i=size(X,2):-1:1  % we want to count backwards here
    % Apply constraints
    isPhysical=applyConstraints(X(:,i),extPar);
    disp('------ Fitness ------');
    if isPhysical
        % Parameters from vector X.
        tic;
        
        % Fixed Parameters
       %{
        ObsTime=extPar.fixed.obstime;
        mass=extPar.fixed.mass;
        radius=extPar.fixed.radius;
        inclination=extPar.fixed.inclination;
        emission=extPar.fixed.emission;
        timeShift=extPar.fixed.phaseshift;
        rho=extPar.fixed.rho;
        temperature=extPar.fixed.spot_temperature;
        distance=extPar.fixed.distance;
        %}
        %disp(mass);
        
        ObsTime=extPar.fixed.obstime;
        
        radius=X(1,i);
        mass=X(2,i);
        inclination=X(3,i);
        emission=X(4,i);
        timeShift=X(5,i);        
        rho=X(6,i);
        temperature=X(7,i);
        distance=X(8,i);
        
        
        %ObsTime=extPar.fixed.obstime;
        %mass=extPar.fixed.mass;
        %radius=extPar.fixed.radius;
        %inclination=extPar.fixed.inclination;
        %emission=extPar.fixed.emission;
        %timeShift=extPar.fixed.phaseshift;
        %rho=extPar.fixed.rho;
        %temperature=extPar.fixed.spot_temperature;
        %distance=extPar.fixed.distance;
        
       
        %display('RRRRadius = ')
        %display(radius)
        
        % Initialize background
        %{
        background1=extPar.fixed.background(1);
        background2=extPar.fixed.background(2);
        background3=extPar.fixed.background(3);
        background4=extPar.fixed.background(4);
        background5=extPar.fixed.background(5);
        background6=extPar.fixed.background(6);
        background7=extPar.fixed.background(7);
        background8=extPar.fixed.background(8);
        background9=extPar.fixed.background(9);
        background10=extPar.fixed.background(10);
        background11=extPar.fixed.background(11);
        background12=extPar.fixed.background(12);
        background13=extPar.fixed.background(13);
        background14=extPar.fixed.background(14);
        background15=extPar.fixed.background(15);
        %}
        
  
        
        for j = 1:300
           par.back(j) = X(8+j,i);
           %disp(par.back(j))
        end
        
        
     
        %disp(background1)
        time=toc

        %NN=extPar.fixed.numbands-1;
        
cmd = '[Fspot(i),auxOutput{i}] = spotMex_new(mass, radius, extPar.fixed.freq, inclination, emission, timeShift, extPar.fixed.numbins, extPar.fixed.modelchoice, rho, temperature, distance, extPar.fixed.numtheta, extPar.fixed.spectral_model, extPar.fixed.numbands, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.fixed.spots_2, extPar.fixed.obsdata.t, extPar.fixed.bend_file_is, extPar.fixed.mr, extPar.fixed.b, extPar.fixed.psi, extPar.fixed.dcosa, extPar.fixed.toa, extPar.fixed.spotshape, ObsTime, extPar.fixed.inst_curve, extPar.fixed.attenuation, extPar.fixed.inte, extPar.fixed.angl, extPar.fixed.energy';
for j = 1:300    
    %cmd = [cmd,', obsdata(:,',num2str(i*2),'), obsdata(:,',num2str(i*2+1),'), background(',num2str(i),')'];
    %cmd = [cmd,', extPar.fixed.obsdata.f(:,',num2str(j),'),extPar.fixed.background(',num2str(j),')'];
     cmd = [cmd,', extPar.fixed.obsdata.f(:,',num2str(j),'),par.back(',num2str(j),')'];
end

for j = 1:400
    cmd = [cmd,', extPar.fixed.instru(',num2str(j),',:)'];
end

cmd = [cmd,');'];
%disp('This is the command');
%disp(cmd)

        
        
%{
        cmd = '[Fspot(i),auxOutput{i}] = spotMex_trial(mass, radius, extPar.fixed.freq, inclination, emission, timeShift, extPar.fixed.numbins, extPar.fixed.modelchoice, rho, temperature, distance, extPar.fixed.numtheta, extPar.fixed.spectral_model, extPar.fixed.numbands, extPar.fixed.E_band_lower_1, extPar.fixed.E_band_upper_1, extPar.fixed.beaming_model, extPar.fixed.spots_2, extPar.obsdata2.t, extPar.fixed.bend_file_is, extPar.fixed.mr, extPar.fixed.b, extPar.fixed.psi, extPar.fixed.dcosa, extPar.fixed.toa, extPar.fixed.spotshape, ObsTime, extPar.fixed.inst_curve, extPar.fixed.attenuation, extPar.fixed.inte, extPar.fixed.angl';
        for j = 1:extPar.fixed.numbands
            cmd = [cmd,', extPar.obsdata2.f(',num2str(j),',:), extPar.obsdata2.err(',num2str(j),',:), background',num2str(j),''];
        end
        cmd = [cmd,');'];
%}
        disp(cmd)
        disp(i);
        %disp(timeShift);
        eval(cmd);
        disp(Fspot(i));
        time=toc
    
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
%t=X(strcmpi(parameterNames,'temperature'));
%m=1.6;
%r=11.8;
%mr_ratio=extPar.fixed.mass/extPar.fixed.radius;
%mr_ratio=m/r;
mr_ratio=1.477*m/r;
%tsquare=t*t;
%tsquare constraint comes from ML2015
%thing=12.0/(5.0*(1.0-2.0*mr_ratio));
%isPhysical=mr_ratio > 0.0677 && mr_ratio < 0.203;
%isPhysical=mr_ratio > 0.1 && mr_ratio < 0.284 && tsquare<thing+0.1 && tsquare > thing-0.1;
isPhysical=mr_ratio > 0.11 && mr_ratio < 0.284;
%isPhysical=mr_ratio > 0.0677 && mr_ratio < 0.175;
