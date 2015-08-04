
% loads hot wire bridge voltages, uses calibrations to return velocities.  
%loadvelocityff.m
% 
% SINGLE CHANNEL VERSION OF LOADVELOCITYNF
% 
% V = loadvelocityff(datafilename, calibrationfilename, probenumber, sensornumber, ...  
%                    TfilenameorT); 
% 
% 
% Greg Bewley  2011

function [E,a,b,n] = loadvelocityff(datafilename, calibrationfilename, probenumber, sensornumber)

tstart = tic; 
 
  % ---- read the raw data: 
fprintf('  loading the raw data from %s...  \n', datafilename); 
E = loadLabviewbin(datafilename); 

probeinfo = extractprobeinfo(calibrationfilename); 

E = E/probeinfo.probe(probenumber).sensor(sensornumber).gain + ...  
    probeinfo.probe(probenumber).sensor(sensornumber).offset;

if nargin > 4  % NOT IMPLEMENTED YET
	  % ---- convert temperature to degrees C, if necessary: 
	fprintf('  converting the temp...  mean, std temp for each block: \n'); 
	tic
	if (T{1}(1) < 5)  % if the value is less than 5, its probably a voltage, not a temp.  
		T{j} = convertDantectemp(T{j}); 
	end
	T0(j)  =  mean(T{j}); 
	Tstd(j) = std(T{j}); 
	fprintf('%.2f C, %.2f C;  ', T0(j), Tstd(j)); 
	fprintf('\n'); 
	
	  % all data will be corrected to the mean temperature for this collection of data: 
	T0 = mean(T0); 
	toc
	
	  % ---- correct voltages for temperature: 
	fprintf('  correcting voltages for temps...  \n'); 
	tic
	calibdata = extractcalibdata(calibrationfilename); 
	E{j}(i, :) = applytempcorrection(E{j}(i, :), mean(T{j}), calibdata.Twire(i), T0); 
	toc
else
	T0 = []; 
end

%disp(mean(E));
  % ---- process the calibration data: 
fprintf('  processing the calibrations...  \n'); 
[a b n] = processcalibrationff(calibrationfilename, T0, true, probenumber, sensornumber);
fprintf('    a = %0.3f, b = %0.3f, n = %0.3f\n', a, b, n); 

  % ---- apply the calibration:
fprintf('  converting voltages to velocities...  \n'); 
%E = applyKingslaw(E, a, b, n);
E = (a*E.^2)+(b*E)+n;
  % E is now velocity.
totaltime = toc(tstart); 
fprintf('  total time to load and convert the data: %.1f seconds.  \n', totaltime); 
