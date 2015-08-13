% Client program for using Edecfit
% This program takes data from multiple .mat workspaces, and extracts
% necessary parameters for the Edecfit.m energy decay fit function
% Written by Nathan Wei, summer 2015
% Requirements: any number of .mat workspaces saved in a folder, with names
% in format '...n.mat' (where n is the test number)
% Dependencies: Edecfit.m

% -------------------------- PARAMETERS TO SET --------------------------
pathname = fileparts('/Users/nathan/Documents/Data/data08_13_15/');
%filebase = 'statscorr_th2.6th2_0813_0'; % files numbered 0 to tests-1
filebase = 'statscorr_th2.6th2_0813_0';
tests = 10; % number of data collection points along the tunnel
% record distance from active grid to probe (in meters) for each data set
% note: the first distance should correspond to test 0, the 2nd to test 1,
% and so on
dist = [1.6356; 1.764634032; 1.948293331; 2.209702958; 2.581777736; ...
    3.111366651; 3.865151737; 4.938044201; 6.465134961; 8.638704324];
% initial parameters for curve fit (suggested: 1 3 -1)
b0 = [1, 3, -1];
% -----------------------------------------------------------------------

addpath(pathname);
tic;
% validation
if (length(dist) ~= tests)
    error('Number of distances does not match number of tests!');
end

fprintf('Loading data from .mat files... \n');

U_avg = zeros(tests, 1); % array for storing mean velocities
eps_avg = zeros(tests, 1); % array for storing mean energy (U^2)

for i = 1 : tests
    filestr = strcat(filebase, num2str(i-1), '.mat'); % find file name
    disp(filestr);
    vars = load(filestr,'MASvsm','MASvss');
    U_avg(i) = vars.MASvsm; % mean flow velocity
    eps_avg(i) = vars.MASvss^2; % rms velocity squared = average energy
end

fprintf(' done in %.1f seconds. \nSending data to Edecfit for processing.', round(10*toc)/10);

result = Edecfit(dist, U_avg, eps_avg, b0);

