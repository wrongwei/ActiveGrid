%{
Load in velocities from data file and check isotropy
Improved version that uses a full velocity calibration (i.e. multiple 
velocity calibrations done at different angles) for x-wire probe

Nathan Wei, Willem van de Water, and Kevin Griffin, 21 August 2015

Dependencies: loadvelocityxwire.m, processcalibrationxwire.m
Requirements: an x-wire calibration file and n velocity calibration .txt
files (saved directly from StreamWare software, one for each angle 1->n)
Outputs std(streamwisevelocity)/std(normalvelocity), which should be = 1 in
the isotropic case
Keep in mind which channel is which wire (don't reverse them!)
%}

% ----------- PARAMETERS TO CHANGE (for standalone operation) -------------
path = '/Users/nathan/Documents/Data/data08_21_15/'; % location of calib file
folder = 'th3.9th3_iso1'; % name of folder containing data
calibfile = 'calib8_21b.m'; % calibration file name
temperature = 22.2; % average temperature of run (can be left as [] if not taken)
% -------------------------------------------------------------------------
begintime = tic;
addpath(fileparts(path));
addpath(fileparts(strcat(path,folder,'/')));

% convert voltages from .dat files to velocities using Willem's method
fprintf('Loading x-wire velocities:\n');
[ux1 uy1 s1] = loadvelocityxwire('xpos100_ypos100_evts0-2999999XCh1_Ch1.dat',...
    'xpos100_ypos100_evts0-2999999XCh2_Ch4.dat', calibfile, 1, temperature);
ux = [ux1];
uy = [uy1];
sorry = s1;
total = length(ux);

ux(isnan(ux)) = []; % remove NaN parameters from ux
uy(isnan(uy)) = []; % same, assuming uy NaNs occur in same places ux ones did

iso = std(ux)/std(uy);
endtime = toc(begintime);
fprintf('Finished loading and converting all data sets in %.1f seconds.\n\n', round(10*endtime)/10); 
fprintf('Mean stream velocity = %f m/s \n', mean(ux));
fprintf('Mean normal velocity = %f m/s \n', mean(uy));
fprintf('Coefficient of isotropy = %f \n', iso);
fprintf('Percentage of problematic data points = %.2f%%\n', (sorry/total)*100);

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