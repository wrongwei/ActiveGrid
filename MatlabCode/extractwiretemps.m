
% extract wire temperatures from a calibration file.  
%extractwiretemps.m
% 
% Twire = extractwiretemps(calibrationfilename); 
% 
% Greg Bewley 2010

function Twire = extractwiretemps(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry);
