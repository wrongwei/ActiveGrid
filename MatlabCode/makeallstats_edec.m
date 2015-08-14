% Script to automatically process data files from an energy decay data
% collection set, currently built for 10 minutes of data at each point
% Uses makeallstats (as a function) multiple times to save workspaces
% Requires: a directory with calibration file, and a number of data folders
%  indexed by test number, located either within the directory or in a
%  sub-folder in the directory
% Dependencies: makeallstats
% Nathan Wei, August 2015

% -------------------------- PARAMETERS TO SET --------------------------
foldername = '/Users/nathan/Documents/Data/data08_13_15/';
tests = 10; % number of data collection points along the tunnel
calibfile = 'calib8_13.m';
% generic name of folder holding each data point - can include subfolder
% here if it exists. Ex. 'folder0' for data folder in directory, or
% 'subfolder/datafolder0' for data folder in sub-folder.
datafolderbase = 'th2.6th2/';
% standardized beginning of each workspace makeallstats will create
outputfilebase = 'statscorr_th2.6th2_0813_0';
% temperatures measured at each point - leave empty if no temperature data was taken
testtemps = [22.62; 22.635; 22.65; 22.655; 22.68; 22.7; 22.71; 22.72; 22.735; 22.75];
%testtemps = [22.74; 22.75; 22.745; 22.73; 22.72; 22.71; 22.71; 22.735; 22.735; 22.735];
needU = true; % save vector u in workspace if you want error bars in EdecfitClient.m
% -----------------------------------------------------------------------

tstart = tic;

for t = 1 : tests
    disp(strcat('Processing test #', num2str(t), '/', num2str(tests)));
    folderstring = strcat(datafolderbase, num2str(t-1),'/');
    outputstring = strcat(outputfilebase, num2str(t-1), '.mat');
    if isempty(testtemps)
        temp = []; % so makeallstats knows not to do temp correction
    else
        temp = testtemps(tests);
    end
    makeallstats(foldername, folderstring, calibfile, outputstring, temp, needU); 
    % note: makeallstats adds the final index and file extension
end

% play sound to alert sleeping user to end of data processing
t = 0:(1/8000):0.25;
y1 = sin(2*pi*440*t);
y2 = sin(2*pi*660*t);
y3 = sin(2*pi*880*t);
y = [y1 y2 y3 y1];
sound(y, 8000);

fprintf('\n\nTotal time to create %d workspaces: %.1f seconds.\n', tests, round(10*toc(tstart))/10);