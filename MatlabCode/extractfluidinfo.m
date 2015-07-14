
% extract fluid info from a calibration file.  
%extractprobeinfo.m
% 
% fluidinfo = extractfluidinfo(calibrationfilename); 
% 
% fluidinfo.medium   % the gas as a string.  
% fluidinfo.p        % [mbar] actual pressure of the gas.  
% fluidinfo.T        % [C] temperature of the fluid.  
% fluidinfo.nu       % [m^2/s] fluid viscosity.  
% 
% Greg Bewley 2011

function fluidinfo = extractfluidinfo(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd; 
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry); 
