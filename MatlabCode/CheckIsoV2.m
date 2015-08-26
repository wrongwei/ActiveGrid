%{
---------------------------------------------------------------------------
Load in velocities from data file and check isotropy
2D version of loadvelocityff that uses a full velocity calibration (i.e. 
several velocity calibrations done at different angles) for an x-wire probe
Calculates the mean streamwise and normal velocities from a given data set,
the standard deviations (rms velocities) in both directions, and the
isotropy coefficient (std(vs)/std(vn))
Also displays the percentage of points that were not able to be resolved by
the calibration algorithm in loadvelocityxwire.m

Nathan Wei, Willem van de Water, and Kevin Griffin
21 August 2015

Dependencies: loadvelocityxwire.m, processcalibrationxwire.m
Requirements: an x-wire calibration file and n velocity calibration .txt
files (saved directly from StreamWare software, one for each angle 1->n)
Outputs std(streamwisevelocity)/std(normalvelocity), which should be = 1 in
the isotropic case
Keep in mind which channel is which wire (don't reverse them!)
---------------------------------------------------------------------------
%}

% ------------------------ PARAMETERS TO CHANGE ---------------------------
path = '/Users/nathan/Documents/Data/data08_21_15/'; % location of calib file
folder = 'lt5.2lt50_isodecay/9'; % name of folder containing data
calibfile = 'calib8_21b.m'; % calibration file name
temperature = 21.93; % average temperature of run (can be left as [] if not taken)
% -------------------------------------------------------------------------
begintime = tic;
addpath(fileparts(path));
addpath(fileparts(strcat(path,folder,'/')));

% convert voltages from .dat files to velocities using Willem's method
fprintf('Loading x-wire velocities:\n');
% run computations on each file pair, then combine resulting vectors
[ux1 uy1 s1] = loadvelocityxwire('xpos100_ypos100_evts0-2999999XCh1_Ch1.dat',...
    'xpos100_ypos100_evts0-2999999XCh2_Ch4.dat', calibfile, 1, temperature);
[ux2 uy2 s2] = loadvelocityxwire('xpos100_ypos100_evts3000000-5999999XCh1_Ch1.dat',...
    'xpos100_ypos100_evts3000000-5999999XCh2_Ch4.dat', calibfile, 1, temperature);
[ux3 uy3 s3] = loadvelocityxwire('xpos100_ypos100_evts6000000-8999999XCh1_Ch1.dat',...
    'xpos100_ypos100_evts6000000-8999999XCh2_Ch4.dat', calibfile, 1, temperature);
[ux4 uy4 s4] = loadvelocityxwire('xpos100_ypos100_evts9000000-11999999XCh1_Ch1.dat',...
    'xpos100_ypos100_evts9000000-11999999XCh2_Ch4.dat', calibfile, 1, temperature);
ux = [ux1;ux2;ux3;ux4];
uy = [uy1;uy2;uy3;uy4];
sorry = s1+s2+s3+s4; % total number of bad points
total = length(ux); % total number of points (before removing bad points)

ux(isnan(ux)) = []; % remove NaN parameters from ux
uy(isnan(uy)) = []; % same, assuming uy NaNs occur in same places ux ones did

stdx = std(ux);
stdy = std(uy);
iso = stdx/stdy; % calculate isotropy coefficient
endtime = toc(begintime);

% print statistics
fprintf('Finished loading and converting all data sets in %.1f seconds.\n\n', round(10*endtime)/10); 
fprintf('Mean stream velocity = %f m/s \n', mean(ux));
fprintf('Mean normal velocity = %f m/s \n', mean(uy));
fprintf('RMS stream velocity = %f m/s \n', stdx);
fprintf('RMS normal velocity = %f m/s \n', stdy);
fprintf('Coefficient of isotropy = %f \n', iso);
fprintf('Percentage of problematic data points = %.2f%%\n', (sorry/total)*100);
%fprintf('%f\t%f\t%f\t%f\t%f\n', mean(ux),mean(uy),stdx,stdy,iso);

% play sound to alert sleeping user to end of data processing
% (I'm having way too much fun with this)
t1 = 0:(1/8000):0.54;
t2 = 0:(1/8000):0.18;
E1 = sin(2*pi*329.63*t2);
p = zeros(1,100);
A1 = sin(2*pi*440*t1);
B = sin(2*pi*493.88*t2);
Cs = sin(2*pi*554.37*t2);
D = sin(2*pi*587.33*t2);
E2 = sin(2*pi*659.25*t1);
A2 = sin(2*pi*880*t1);
y = [E1 p E1 p E1 p A1 A1 E2 E2 D Cs B A2 A2];
sound(y, 8000);

rmpath(fileparts(path));
rmpath(fileparts(strcat(path,folder,'/')));