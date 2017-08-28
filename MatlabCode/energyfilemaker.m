% This program is a misnomer. It differs from EdecfitClient in that it does
% not call Edecfit
% This program takes data from multiple statscorr.mat workspaces, and extracts
% necessary parameters to make an energy vs  distance along tunnel graph 
% Written by Kevin Griffin and Nathan Wei, summer 2015
% Requirements: any number of .mat workspaces saved in a folder, with names
%  in format '...n.mat' (where n is the test number). NOTE: make sure you
%  comment out the line in makeallstats that clears u from the workspace!
% Dependencies: Edecfit.m

% -------------------------- PARAMETERS TO SET --------------------------
pathname = fileparts('/Users/kevin/Documents/Data/data08_27_15/');
%filebase = 'statscorr_th2.6th2_0813_0'; % files numbered 0 to tests-1
filebase = '/statscorr_lt5.2lt4_0827_0';
energyFileStr = '/energy_lt5.2lt4_0827';
tests = 10; % number of data collection points along the tunnel
startingTestNumberPlus1 = 1;
%startingTestNumberPlus1 = 6;
% record distance from active grid to probe (in meters) for each data set
% note: the first distance should correspond to test 0, the 2nd to test 1,
% and so on

%dist for 0813 spacing
%dist = [9.638704324;7.465134961;5.938044201;4.865151737;4.111366651;3.581777736;3.209702958;2.948293331;2.764634032;2.6356];
%dist for 0814 spacingf
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

%U_avg = zeros(tests, 1); % array for storing mean velocities
%eps_avg = zeros(tests, 1); % array for storing mean energy (U^2)

%stderr = zeros(tests, 1); % array for calculating and storing standard error
energyArray = zeros(tests, 1);
%distanceArray = zeros(tests, 1);
stdErrArray = zeros(tests, 1);
energyVarianceArray = zeros(tests, 1);
oneOverEScaleArray = zeros(tests, 1); % array for storing corr. function length scales
MASvsmArray = zeros(tests,1);

for i = startingTestNumberPlus1 : tests
    filestr = strcat(pathname,filebase, num2str(i-1), '.mat'); % find file name
    fprintf('Currently processing %s\n',filestr);
    tic
    if (errorbars)
        fprintf('Loading Data...');
        vars = load(filestr,'MASvsm','MASvss','rCmax','oneOverEScale','u');
        % error of velocity flucuations
        %stderr(i) = vars.MASvss*sqrt(vars.oneOverEScale*samplingFrequency/vars.rCmax/vars.MASvsm);
        %error of energy normalized by mean velocity squared        
        numberOfCorrLengths = vars.rCmax * vars.MASvsm / vars.oneOverEScale / samplingFrequency;
        energyVarianceArray(i) = mean(((vars.u-vars.MASvsm).^2 - vars.MASvss^2).^2);
        stdErrArray(i) = sqrt(energyVarianceArray(i)/numberOfCorrLengths);
        %normalize
        stdErrArray(i) = stdErrArray(i)/vars.MASvsm^2;
    else
        fprintf('Loading Data...');
        vars = load(filestr,'MASvsm','MASvss','rCmax','oneOverEScale');
    end
    fprintf(' Done in %.1f seconds.\nPerforming energy calculations...', round(10*toc)/10);
    tic
    oneOverEScaleArray(i) = vars.oneOverEScale;
    paddled = 0.115; % distance between two adjacent paddles (m)
    %normalize
    energyArray(i) = vars.MASvss^2;
    MASvsmArray(i) = vars.MASvsm;
    fprintf('Done in %.1f seconds.\nSaving Data...', round(10*toc)/10);
end
fprintf('Plotting...');
%figure
if (errorbars)
    errorbar(dist/paddled,energyArray./MASvsmArray.^2, stdErrArray,'.');
else
    scatter(dist/paddled,energyArray./MASvsmArray.^2, 1000, '.');
end
hold on;

tic
fprintf('Saving...');
energyFileStr = strcat(pathname,energyFileStr,'.mat');
rCmax = vars.rCmax;
if(errorbars)
    save(energyFileStr,'dist','energyArray','rCmax','oneOverEScaleArray','MASvsmArray','energyVarianceArray');
else
    save(energyFileStr,'dist','energyArray','rCmax','oneOverEScaleArray','MASvsmArray');
end
fprintf('Done in %.1f seconds.\n', round(10*toc)/10);


% plot 1/e scale against distance down tunnel
figure;
plot(dist,oneOverEScaleArray);
xlabel('Normalized Distance');
ylabel('1/e Length Scale');
title('Evolution of 1/e Length Scale with Distance along Tunnel');