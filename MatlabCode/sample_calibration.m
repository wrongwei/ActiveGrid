
% experiment hot wire calibration file
% 
% output is at least a vector with fluid velocities, 
% and a vector with bridge voltages.  

clear exptinfo fluidinfo probeinfo calibdata

  % ---- information about the setup: 
  % info about the experiment: 
exptinfo.date   = [2015 8 5];  % date of the experiment [year month day].  
exptinfo.medium = 'Air';         % the gas as a string.  
exptinfo.pnom   = 1000;          % [mbar] nominal pressure of the gas.  

  % ---- information about the fluid: 
fluidinfo.medium = 'Air';        % the gas as a string.  
fluidinfo.p      = 1000;         % [mbar] actual pressure of the gas.  
fluidinfo.T      = nan;          % [C] temperature of the fluid.  
dens = 1.2041; 
visc = 15.11*10^-6; 
fluidinfo.nu     = visc;         % [m^2/s] fluid viscosity.  
fluidinfo.dens   = dens;         % [kg/m^3] fluid density.  

  % ---- info about the probes: 
probeinfo.notes = '1 P11 wire.  '; 
probeinfo.samplerate = 20;       % [kHz] rate at which samples were collected.  
probeinfo.numsmpls = 512;      % number of samples collected at each calibration point.  

  % info about each probe: 
probeinfo.probe(1).probename = 'S'; 

  % info about each sensor: 
probeinfo.probe(1).sensor(1).channel    = 1; 
probeinfo.probe(1).sensor(1).resistance = 3.516;  % [Ohms]
probeinfo.probe(1).sensor(1).offset     = 1.577;  % [Volts]
probeinfo.probe(1).sensor(1).gain       = 64;
probeinfo.probe(1).sensor(1).overheat   = 0.8;
% temperature coefficient of resistance for wire
probeinfo.probe(1).sensor(1).alpha      = 3.6e-3;

% velocity calibration .txt file (6 columns, we're interested in 1+2)
fileID = fopen('velocity_calibration_8_05.txt');
data = textscan(fileID,'%f %f %f %f %f %f');
fclose(fileID);

% UNCOMMENT THIS LINE if you want to use a polynomial fit
calibdata.polyfit = 2; % order of fit

  % ---- derived quantities: 
  % number of probes: 
probeinfo.numprobes = length(probeinfo.probe); 
  % number of sensors on each probe: 
for icalib = 1:probeinfo.numprobes
	probeinfo.probe(icalib).numsensors = length(probeinfo.probe(icalib).sensor); 
end
  % total number of channels used: 
n = 0; for icalib = 1:probeinfo.numprobes, n = n + probeinfo.probe(icalib).numsensors; end
probeinfo.numchannels = n; 

  % ---- angular calibration: 
  % probe wire angles, in degrees: 
probeinfo.probe(1).sensor(1).angle = 0; 


  % ---- temperature calibration: 
  % by CHANNEL: 
calibdata.Twire = [nan]; 

  % copy effective wire temperature in degrees C for each wire: 
for icalib = 1:probeinfo.numprobes
	for jcalib = 1:probeinfo.probe(icalib).numsensors
		probeinfo.probe(icalib).sensor(jcalib).Twire = ...  
		    calibdata.Twire(probeinfo.probe(icalib).sensor(jcalib).channel); 
	end
end


  % ---- velocity calibration, all values in order of CHANNEL: 
  % fan speed in Hz, or in percent: 
calibdata.fanspeed  =  []; 

  % Pitot tube pressure difference [mbar]: 
poffset = 0.00;   % [mbar] *extrapolated* zero value of differential pressure transducer
calibdata.dp     =     []; % - poffset; 

  % fluid density [kg/m^3]: 
calibdata.dens    =    [fluidinfo.dens]; 

  % ---  NSTAP  ---
  % fluid temperature [C]: 
calibdata.T(1, :)  =  0.5*(21.8 + 22.9); 

% fluid velocity [m/s]: (if you want Matlab to do curve fitting)
calibdata.v(1, :)  =  cell2mat(data(1));
% calibration potential measurements [V]
calibdata.E(1, :)  =  cell2mat(data(2));
% calibration temperatures [C]
calibdata.calibtemp(1, :) = cell2mat(data(3));

% Otherwise, use synthesized calibration curves from Dantec coefficients
%fprintf('warning: synthesized calibration'); 
%{
% !!!  synthesized calibration data from Dantec power law fit coefficients !!!  
calibdata.v(1, :)  =  [1:0.1:2]; % velocity range

% wire voltages [V]: (power law coefficients from Dantec)
Abef = -52.373753;  Bbef = 52.051975;  nbef = 1.01; 
Aaft = -52.373753;  Baft = 52.051975;  naft = 1.01; 

% polynomial fit coefficients (from Excel or Dantec)
C2 = -3.2697; C1 = 15.815; C0 = -11.265;

% generate voltage curve from velocity range specified above
%calibdata.E(1, :)  =  0.5*(sqrt(Abef + Bbef*calibdata.v(1, :).^nbef) + ...  
%                           sqrt(Aaft + Baft*calibdata.v(1, :).^naft)); 
calibdata.E(1, :)  = (C2*calibdata.v(1, :).^2)+(C1*calibdata.v(1, :))+C0;
calibdata.E(1, :)  = calibdata.E(1, :); 
%}