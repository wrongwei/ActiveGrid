
% to extract the hot wire probe angles from a calibration file.  
%extractangles.m
% 
% probeangles = extractangles(calibrationfilename); 
% 
% 
% Greg Bewley 2010

function probeangles = extractangles(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry);
