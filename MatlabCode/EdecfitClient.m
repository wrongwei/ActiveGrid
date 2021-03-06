%{
---------------------------------------------------------------------------
Client program for using Edecfit
This program takes data from multiple .mat workspaces and extracts
necessary parameters for the Edecfit.m energy decay fit function.
This also calculates and stores error bar magnitudes (calculated from the
standard error of the normalized energy)

Written by Nathan Wei and Kevin Griffin, summer 2015

Requirements: any number of .mat workspaces saved in a folder, with names
 in format '...n.mat' (where n is the test number). NOTE: make sure you
 comment out the line in makeallstats that clears u from the workspace!
Dependencies: Edecfit.m
---------------------------------------------------------------------------
%}

% -------------------------- PARAMETERS TO SET --------------------------
pathname = fileparts('/Users/nathan/Documents/Data/data08_19_15/');
%filebase = 'statscorr_th2.6th2_0813_0'; % files numbered 0 to tests-1
filebase = 'statscorr_lt5.2lt50_0819_0';
tests = 10; % number of data collection points along the tunnel
startingTestNumberPlus1 = 1; % test to start analysis at, plus 1
%startingTestNumberPlus1 = 6;
% record distance from active grid to probe (in meters) for each data set
% note: the first distance should correspond to test 0, the 2nd to test 1,
% and so on

%dist for 0813 spacing
%dist = [9.638704324;7.465134961;5.938044201;4.865151737;4.111366651;3.581777736;3.209702958;2.948293331;2.764634032;2.6356];
%dist for 0814 spacing (in meters from grid)
dist = [9.638704324;8.345409279;7.225645035;6.256127701;5.416697557;4.68989986;4.060621894;3.515778729;3.044041133;2.6356];
% dist for tests 0-4
%dist = [9.638704324;8.345409279;7.225645035;6.256127701;5.416697557];
% dist for tests 5-9
%dist = [4.68989986;4.060621894;3.515778729;3.044041133;2.6356];

% initial parameters for curve fit (suggested: 1 3 -1)
b0 = [1, 3, -1];
errorbars = true; % calculate and plot error bars on graphs if true
samplingFrequency = 20000; % hz
paddled = 0.115; % distance between two adjacent paddles (m)
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
oneOverEscales = zeros(tests, 1); % array for storing corr. function length scales
stderr = zeros(tests, 1); % array for calculating and storing standard error

for i = startingTestNumberPlus1 : tests
    filestr = strcat(filebase, num2str(i-1), '.mat'); % find file name
    disp(filestr);
    vars = load(filestr,'MASvsm','MASvss','rCmax','oneOverEScale','u');
    % error of energy normalized by mean velocity squared        
    if (errorbars)
        numberOfCorrLengths = vars.rCmax * vars.MASvsm / vars.oneOverEScale...
            / samplingFrequency; % effective number of samples (corr. lengths)
        % standard deviation of energy (velocity fluctuations squared,
        % compared to the RMS fluctuation velocity)
        stddev = sqrt(mean(((vars.u-vars.MASvsm).^2 - vars.MASvss^2).^2));
        % calculate standard error using above parameters
        stderr(i) = stddev / sqrt(numberOfCorrLengths);
        stderr(i) = stderr(i)/vars.MASvsm^2 % normalize by mean velocity ^2
    end
    U_avg(i) = vars.MASvsm; % mean flow velocity
    eps_avg(i) = vars.MASvss^2; % rms velocity squared = average energy
    oneOverEscales(i) = vars.oneOverEScale; % 1/e scale
end

% call Edecfit.m with these calculated parameters
fprintf(' done in %.1f seconds. \nSending data to Edecfit for processing.', round(10*toc)/10);
if (errorbars)
    result = Edecfit(dist, U_avg, eps_avg, b0, stderr);
else
    result = Edecfit(dist, U_avg, eps_avg, b0);
end

% plot 1/e scale against distance down tunnel
figure;
plot(dist,oneOverEscales);
xlabel('Normalized Distance');
ylabel('1/e Length Scale');
title('Evolution of 1/e Length Scale with Distance along Tunnel');
