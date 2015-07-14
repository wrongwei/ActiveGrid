
% load new hot wire data format
%loadLabviewbin.m
% 
% E = loadLabviewbin(filename); 
% 
% each file contains a list of voltages from one channel.  
% 
% Greg Bewley  Aug. 2011

function E = loadLabviewbin(filename)

offset = 5; 
%amplification = 6353;  % error in original Labview code.  
amplification = 6553; 

  % open the file: 
fid = fopen(filename, 'r', 'ieee-le'); 

if (fid == -1)
	fprintf('loadLabviewbin: problem opening the file %s.  \n', filename); 
	E = nan; 
	
else
	  % the data are signed 16 bit long endian integers: 
	E = fread(fid, inf, 'int16'); 
	E = E/amplification + offset; 
	
	fclose(fid); 
end
