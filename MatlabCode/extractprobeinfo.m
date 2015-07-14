
% extract probe info from a calibration file.  
%extractprobeinfo.m
% 
% probeinfo = extractprobeinfo(calibrationfilename); 
% 
% probeinfo: 
% .notes
% .samplerate
% .numsmpls
% .numchannels
% .numprobes
% .probe
% 
% probe: 
% .probename
% .numsensors
% .sensor
% 
% sensor: 
% .channel
% .angle
% .Twire
% 
% Greg Bewley 2010

function probeinfo = extractprobeinfo(calibrationfilename)

[pathstr, filename] = fileparts(calibrationfilename); 

odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(filename); 
cd(odiry);
