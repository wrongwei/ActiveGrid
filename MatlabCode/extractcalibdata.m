
% extract wire temperatures from a calibration file.  
%extractcalibdata.m
% 
% calibdata = extractcalibdata(calibrationfilename); 
% 
% calibdata: 
% .Twire
% .fanspeed
% .dp
% .dens
% .T  -  fluid temperature
% .v  -  fluid velocity
% .E  -  wire voltage
% 
% Greg Bewley 2010

function calibdata = extractcalibdata(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry);
