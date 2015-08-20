% Script to take .mat files from multiple iterations of a data test & merge
% them into one giant .mat file (carrying out makeallstats.m calculations)
% Requirements: .mat files including the parameters 'u' and 'oneOverEScale'
% in the directory specified by 'path'
% Nathan Wei, August 2015
clear all;

% -------------------------- PARAMETERS TO SET -------------------------- %
path = '/Users/nathan/Documents/Data/lt5.2lt50_all/'; % path
searchstring = '*09.mat'; % regular expression for files you want to merge
outputname = 'lt5.2lt50_merged_09.mat';
% ----------------------------------------------------------------------- %
addpath(fileparts(path));

fprintf('Extracting and combining velocities from these files:\n');
tic;

files = dir(strcat(path, searchstring)); % find all files matching searchstring
oneOverEScale = 0;

%extract velocity and average 1/e scale
u = [];
for i = 1 : length(files)
    disp(files(i).name); % debugging
    temp = load(strcat(path, files(i).name),'u','oneOverEScale');
    oneOverEScale = oneOverEScale + (temp.oneOverEScale*length(temp.u)); % weighted average
    u = cat(1, u, temp.u);
end
clear files; % don't let these get into our workspace
clear temp;
fprintf(' Velocities merged in %.1f seconds.\nRunning computations...', round(10*toc)/10); 

% double check that u is a single and not a double
u = single(u);

% number of correlation function separations starting from one in samples: 
rCmax = single(length(u));

tic
% compute mean: (weird names are carry-over from makeallstats.m)
MASvsm = single(mean(u)); 

% compute standard deviation / RMSD: 
MASvss = single(std(u));

% compute 1/e scale from weighted average of previous values
oneOverEScale = oneOverEScale / rCmax;

fprintf(' done in %.1f seconds.\nSaving data...', round(10*toc)/10); 
tic;

matfile = fullfile(fileparts(path), outputname);
save(matfile);

rmpath(fileparts(path));

% play sound to alert sleeping user to end of data processing
t = 0:(1/8000):0.075;
y1 = sin(2*pi*440*t);
y2 = sin(2*pi*554.37*t);
y3 = sin(2*pi*659.25*t);
y4 = sin(2*pi*880*t);
yfull = [y1 y2 y3 y4 y3 y2 y1 y2 y3 y4 y3 y2 y1];
sound(yfull, 8000);

fprintf(' done in %.1f seconds.\n', round(10*toc)/10);
%close all; 

