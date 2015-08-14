% Client program for using Edecfit
% This program takes data from multiple .mat workspaces, and extracts
% necessary parameters for the Edecfit.m energy decay fit function
% Written by Nathan Wei, summer 2015
% Requirements: any number of .mat workspaces saved in a folder, with names
%  in format '...n.mat' (where n is the test number). NOTE: make sure you
%  comment out the line in makeallstats that clears u from the workspace!
% Dependencies: Edecfit.m

% -------------------------- PARAMETERS TO SET --------------------------
pathname = fileparts('/Users/nathan/Documents/Data/data08_13_15/');
%filebase = 'statscorr_th2.6th2_0813_0'; % files numbered 0 to tests-1
filebase = 'statscorr_th2.6th2_0813_0';
tests = 10; % number of data collection points along the tunnel
% record distance from active grid to probe (in meters) for each data set
% note: the first distance should correspond to test 0, the 2nd to test 1,
% and so on
dist = [8.638704324; 6.465134961; 4.938044201; 3.865151737; 3.111366651;...
    2.581777736; 2.209702958; 1.948293331; 1.764634032; 1.6356];
%dist = [8.638704324; 8.638704324; 6.465134961; 6.465134961; 4.938044201; 4.938044201; 3.865151737; 3.865151737; 3.111366651; 3.111366651;...
%    2.581777736; 2.581777736; 2.209702958; 2.209702958; 1.948293331; 1.948293331; 1.764634032; 1.764634032; 1.6356; 1.6356];
% initial parameters for curve fit (suggested: 1 3 -1)
b0 = [1, 3, -1];
errorbars = true; % calculate and plot error bars on graphs if true
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
stderr = zeros(tests, 1); % array for calculating and storing standard error

for i = 1 : tests
    filestr = strcat(filebase, num2str(i-1), '.mat'); % find file name
    disp(filestr);
    if (~errorbars)
        vars = load(filestr,'MASvsm','MASvss');
    else
        vars = load(filestr,'MASvsm','MASvss','u');
        for x = 1 : length(vars.u) % variance of square of velocity fluctuations
            stderr(i) = stderr(i) + (((vars.u(x)-vars.MASvsm)^2-vars.MASvss^2)^2/length(vars.u));
        end
        stderr(i) = sqrt(stderr(i)); % standard deviation of velocity fluctuations squared
        stderr(i) = stderr(i)/sqrt(length(vars.u)); % standard error of U^2
        stderr(i) = stderr(i)/(vars.MASvsm^2); % error normalized by mean velocity squared
        %disp(stderr(i));
    end
    U_avg(i) = vars.MASvsm; % mean flow velocity
    eps_avg(i) = vars.MASvss^2; % rms velocity squared = average energy
end

fprintf(' done in %.1f seconds. \nSending data to Edecfit for processing.', round(10*toc)/10);
if (errorbars)
    result = Edecfit(dist, U_avg, eps_avg, b0, stderr);
else
    result = Edecfit(dist, U_avg, eps_avg, b0);
end
