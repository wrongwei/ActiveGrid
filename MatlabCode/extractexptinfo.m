
% extract wire temperatures from a calibration file.  
%extractexptinfo.m
% 
% exptinfo = extractexptinfo(calibrationfilename); 
% 
% exptinfo: 
% .date
% .medium
% .pnom
% 
% Greg Bewley 2011

function exptinfo = extractexptinfo(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry);
