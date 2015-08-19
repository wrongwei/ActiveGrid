
% script to return standard statistics of a hot wire signal
%makeallstats.m
% can be run as script or called as function (implemented summer 2015)
% 
% must exis% u      - the velocity
% deltaT - the intersample time
% rCmax  - the upper limit in samples of separations for the correlation function
% rs     - the separations for the structure functions
% histX  - the bin edges for the histogram
% 
% returns: 
% MASvsm - the mean velocity
% MASvss - the standard deviation
% MASp   - the spectrum
% MASf   - the frequencies for MASp
% MASC   - the correlation function
% MASS2, MASS3, MASS3a, MASS4, MASS6 - the structure functions
% MASN   - the histogram
% 

function [] = makeallstats(path, folder, calibfile, outputname, actualtemp, needU)

if (~exist('computestructure')), computestructure = 1; end
if (~exist('highorder')), highorder = 1; end
%fprintf('  working on %d samples.  \n');


% ----------- PARAMETERS TO CHANGE (for standalone operation) -------------
if (nargin == 0) % set up parameters if they're not provided
    %pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/726G0.54/');
    path = '/Users/kevin/Documents/Data/data08_12_15/'; % location of calib file
    folder = 'th2.6th2'; % name of folder containing data
    outputname = 'statscorr_th2.6th2_0812.mat'; % name your .mat workspace!
    calibfile = 'calib8_12.m'; % calibration file name (set here for convenience)
    actualtemp = 22.275; % change this if you have a temperature measurement you want to use, otherwise should be []
    needU = true; % save vector u in workspace - 'false' to make smaller file, 'true' if you need access later
end
% -------------------------------------------------------------------------
samplingFrequency = 20000; % Hz

% add necessary paths and load data files from folder
files = dir(strcat(path, folder, '/*.dat'));
addpath(fileparts(path));
addpath(fileparts(strcat(path, folder, '/')));

%extract velocity
u = [];
for i = 1 : length(files)
    disp(files(i).name); % debugging
    newU = loadvelocityff(files(i).name, calibfile, 1, 1, actualtemp);
    u = cat(1, u, newU);
end
clear newU;
fprintf('velocity extracted.\n');
fprintf('Basic velocity computations... ');

%this is 1/ the sampling frequency
deltaT = 1/samplingFrequency;

%histX = 35;

% number of correlation function separations starting from one in samples: 
rCmax = length(u);

sepval = [1:rCmax]/(samplingFrequency)*mean(u);
 
  % structure function separations in samples: 
%rs = makelogtime(1, 10000, 50); 
%rs = unique(round(rs)); 

tic
  % compute mean: 
MASvsm = mean(u); 

  % compute rms velocity: THIS IS COMPLETELY NUTS - rms = std
%rmsvelocity = rms(u);

  % compute standard deviation / RMSD: 
MASvss = std(u);
fprintf('  done in %.1f seconds.    \n', round(10*toc)/10); 

%figure;
%histogram(u);

  % compute spectrum: 
%tic
%[MASp MASf] = hwspectrum(u-MASvsm, deltaT); 
%fprintf('spectrum...  done in %.1f seconds.\n', round(10*toc)/10); 

  % compute correlation functions (as a check, since already computed
  % through S2):
tic;
fprintf('correlation functions... ');
temp = xcorr(u-MASvsm, rCmax, 'coeff'); 
MASC = single(temp(rCmax + (1:rCmax))); 
fprintf(' done in %.1f seconds.\n', round(10*toc)/10); 

tic

% Find the distance where MASC reaches the correlation value of 1/e.
% We call this distance the oneOverEScale
% We wait till a value less than 1/e has happened twice, just incase
% there is 1 crazy point
fprintf('computing 1/e scale... ');
count = 0; 
n = 2;
for i=1:length(MASC)
    if(MASC(i) <= exp(-1) ) 
        count = count + 1;
    end;
    if(count >= n)
        oneOverEScale = sepval(i); 
        break;
    end;
end;
fprintf(' done in %.1f seconds.\n', round(10*toc)/10);

% compute structure functions: 

tic
%{
if computestructure
	MASS2 = structfunc1mex(u', rs, 2); 
	MASS3a = structfuncabsmex(u', rs, 3); 
	if highorder
		MASS3 = structfunc1mex(u', rs, 3); 
		MASS4 = structfunc1mex(u', rs, 4); 
		MASS6 = structfunc1mex(u', rs, 6); 
	else
		MASS3 = 0; 
		MASS4 = 0; 
		MASS6 = 0; 
	end
else 

	MASS2 = 0; MASS3 = 0; MASS3a = 0; MASS4 = 0; MASS6 = 0; 
end
fprintf('structure functions...  done in %.1f seconds.  histogram...  \n', round(10*toc)/10); 

  % compute histogram: 
tic
MASN = hist(u-MASvsm, histX); 
fprintf('  done in %.1f seconds.  \n', round(10*toc)/10); 

%plot structure functions
H1 = figure(1);
subplot(2,2,1);
loglog(rs,MASS2,'LineWidth', 2);
legend('MASS2','Location','SouthEast');
xlabel('rs');
ylabel('structure function');
subplot(2,2,2);
loglog(rs,MASS3,rs,MASS3a,'LineWidth', 2);
legend('MASS3','MASS3a','Location','SouthEast');
xlabel('rs');
ylabel('structure function');
subplot(2,2,3);
loglog(rs,MASS4,'LineWidth', 2);
legend('MASS4','Location','SouthEast');
xlabel('rs');
ylabel('structure function');
subplot(2,2,4);
loglog(rs,MASS6,'LineWidth', 2);
legend('MASS6','Location','SouthEast');
xlabel('rs');
ylabel('structure function');

%make histogram
H2 = figure(2);
hist(u-MASvsm, histX);

%make correlation function plot 
H3 = figure(3);
scatter(sepval,MASC,'o');
%loglog correlation function plot
H4 = figure(4);
loglog(sepval,MASC,'o');

%save

%}
fprintf('  done in %.1f seconds.\n  Saving data...', round(10*toc)/10); 
tic;
if (~needU)
    clear u; % u is needed in workspace only for error bars on edec runs
end
matfile = fullfile(fileparts(path), outputname);

%structfile = fullfile(pathname, 'struct.fig');
%histfile = fullfile(pathname, 'hist.fig');
%corelfile = fullfile(pathname, 'corel.fig');
%logcorf = fullfile(pathname, 'lcorel.fig');
save(matfile);
%saveas(H1, structfile);
%saveas(H2, histfile);
%saveas(H3, corelfile);
%saveas(H4, lcorelfile);

%save data to three .txt files 


%save('vsm.txt','MASvsm','-append','-ascii');
%save('std.txt','MASvss','-append','-ascii');
%save('rms.txt','rmsvelocity','-append','-ascii');
rmpath(fileparts(path));
rmpath(fileparts(strcat(path, folder, '/')));

% play sound to alert sleeping user to end of data processing
t = 0:(1/8000):0.5;
y = sin(2*pi*440*t);
sound(y, 8000);

fprintf('  done in %.1f seconds.\n', round(10*toc)/10);
%close all; 
