% Client program for using Edecfit
% This program takes data from multiple .mat workspaces, and extracts
% necessary parameters for the Edecfit.m energy decay fit function
% Written by Nathan Wei, 29 July 2015
% Dependencies: Edecfit.m
clear all;
close all;
clc;

% Where is the path?-----------------------------------------------------
pathname = fileparts('/Users/nathan/Documents/Data/data28_07_15/');
addpath(pathname); 
tic;
fprintf('Loading data from .mat files...');
% which workspace stats do you want to load?
% -----------------------------------------------------------------------
files(1) = load('statscorr_28-07-15_g2corr_04ft.mat','MASvsm','u');
files(2) = load('statscorr_28-07-15_g2corr_12ft.mat','MASvsm','u');
files(3) = load('statscorr_28-07-15_g2corr_16ft.mat','MASvsm','u');
files(4) = load('statscorr_28-07-15_g2corr_18ft.mat','MASvsm','u');
files(5) = load('statscorr_28-07-15_g2corr_19ft.mat','MASvsm','u');
files(6) = load('statscorr_28-07-15_g2corr_20ft.mat','MASvsm','u');
% -----------------------------------------------------------------------

% record distances (in meters) from active grid for each data set
% -----------------------------------------------------------------------
dist(1) = 9.4 - 6.854;
dist(2) = 9.4 - 6.5492;
dist(3) = 9.4 - 6.2444;
dist(4) = 9.4 - 5.6348;
dist(5) = 9.4 - 4.4156;
dist(6) = 9.4 - 1.9772;
% -----------------------------------------------------------------------

% enter initial parameters for curve fit (suggested: 1 3 -1)
% -----------------------------------------------------------------------
b0 = [1, 3, -1];
% -----------------------------------------------------------------------
% end of user input requirements

if (length(dist) ~= length(files))
    error('Number of distances must be the same as number of files!');
end
fprintf(' done in %.1f seconds. \nComputing variances...', round(10*toc)/10);
tic;
variance = zeros(1,length(files));
U_avg = zeros(1,length(files));

for i = 1:length(files)
    % create average velocity vector
    U_avg(i) = files(i).MASvsm;
    % calculate variance of velocity over each data set
    for j = 1:length(files(i).u)
        variance(i) = variance(i) + (files(i).u(j) - U_avg(i))^2;
    end
    variance(i) = variance(i) / length(files(i).u); % take average over time
end

fprintf(' done in %.1f seconds. \nSending data to Edecfit for processing.', round(10*toc)/10);

result = Edecfit(dist, U_avg, variance, b0);

