% Script to automatically process data files from an energy decay data
% collection set, currently built for 10 minutes of data at each point
% Uses makeallstats (as a function) multiple times to save workspaces
% Requires: data files named with 2 indices for test number and data set
% number, all in the same folder. A calibration file should be in the
% directory as well.
% Dependencies: makeallstats
% Nathan Wei, August 2015
%{
Example data file name: g2g1_0730_05_3.dat
- g2g1 signifies Gaussians in the spatial and temporal dimensions
- spatial sigma is 2, temporal sigma is 1
- 0730 is the date of the test
- 05 is the position number along the tunnel (0-9 if tests = 10)
- 3 is the data set number (1-4 for 10 minute long data sets)
%}

% -------------------------- PARAMETERS TO SET --------------------------
filename = '/Users/kevin/Documents/Data/data07_31_15/';
tests = 10; % number of data collection points along the tunnel
calibfile = 'calib7_31.m';
datafilebase = 'g3g3_0731_0'; % standardized beginning of each data file
outputfilebase = 'statscorr_g3g3_0731_0'; % standardized beginning of each workspace makeallstats will create
% -----------------------------------------------------------------------

path = fileparts(filename);
addpath(path);
tstart = tic;

for t = 1 : tests
    disp(strcat('Processing test #', num2str(t), '/', num2str(tests)));
    filestring = strcat(datafilebase, num2str(t-1), '_');
    outputstring = strcat(outputfilebase, num2str(t-1));
    makeallstats(filename, calibfile, filestring, outputstring); % note: makeallstats adds the final index and file extension
end

rmpath(path);

fprintf('\n\nTotal time to create %d workspaces: %.1f seconds.\n', tests, round(10*toc(tstart))/10);