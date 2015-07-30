% Client program for using Edecfit
% This program takes data from multiple .mat workspaces, and extracts
% necessary parameters for the Edecfit.m energy decay fit function
% Written by Nathan Wei, summer 2015
% Dependencies: Edecfit.m

% CAN WE JUST USE THE AVERAGE VOLTAGES GIVEN BY THE TXT FILES???

clear all;
close all;
clc;

% Where is the path?-----------------------------------------------------
pathname = fileparts('/Users/nathan/Documents/Data/data07_29_15/');
addpath(pathname); 
tic;
fprintf('Loading data from .mat files...');
% which workspace stats do you want to load?
% -----------------------------------------------------------------------
files(1) = load('statscorr_g2g1_0729_00.mat','MASvsm','u');
files(2) = load('statscorr_g2g1_0729_01.mat','MASvsm','u');
files(3) = load('statscorr_g2g1_0729_02.mat','MASvsm','u');
files(4) = load('statscorr_g2g1_0729_03.mat','MASvsm','u');
files(5) = load('statscorr_g2g1_0729_04.mat','MASvsm','u');
files(6) = load('statscorr_g2g1_0729_05.mat','MASvsm','u');
files(7) = load('statscorr_g2g1_0729_06.mat','MASvsm','u');
files(8) = load('statscorr_g2g1_0729_07.mat','MASvsm','u');
files(9) = load('statscorr_g2g1_0729_08.mat','MASvsm','u');
files(10) = load('statscorr_g2g1_0729_09.mat','MASvsm','u');
% -----------------------------------------------------------------------

% record distances from active grid (in meters) for each data set
% -----------------------------------------------------------------------
dist(1) = 1.6356;
dist(2) = 1.764634032;
dist(3) = 1.948293331;
dist(4) = 2.209702958;
dist(5) = 2.581777736;
dist(6) = 3.111366651;
dist(7) = 3.865151737;
dist(8) = 4.938044201;
dist(9) = 6.465134961;
dist(10) = 8.638704324;
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

