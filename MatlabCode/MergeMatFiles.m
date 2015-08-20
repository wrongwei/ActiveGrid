% Script to take .mat files multiple iterations of a data test and merge
% them into one giant .mat file (carrying out makeallstats.m calculations)
% Nathan Wei, August 2015

% -------------------------- PARAMETERS TO SET -------------------------- %
path = '/Users/nathan/Documents/Data/lt5.2lt50_all/'; % path
searchstring = '*00.mat'; % regular expression for files you want to merge
outputname = 'lt5.2lt50_merged_00.mat';
samplingFrequency = 20000; % Hz
% ----------------------------------------------------------------------- %
addpath(fileparts(path));

fprintf('Extracting and combining velocities from these files:\n');

files = dir(strcat(path, searchstring)); % find all files matching searchstring

%extract velocity
u = [];
for i = 1 : length(files)
    disp(files(i).name); % debugging
    temp = load(strcat(path, files(i).name),'u');
    u = cat(1, u, temp.u);
end
clear files; % don't let these get into our workspace
clear temp;
fprintf('Velocities merged. Basic velocity computations... ');

% ---------- the rest is copied from makeallstats.m ---------- %

% double check that u is a single and not a double
u = single(u);

%this is 1/ the sampling frequency
deltaT = 1/samplingFrequency;

%histX = 35;

% number of correlation function separations starting from one in samples: 
rCmax = single(length(u));

sepval = single([1:rCmax]/(samplingFrequency)*mean(u));
 
  % structure function separations in samples: 
%rs = makelogtime(1, 10000, 50); 
%rs = unique(round(rs)); 

tic
  % compute mean: 
MASvsm = single(mean(u)); 

  % compute rms velocity: THIS IS COMPLETELY NUTS - rms = std
%rmsvelocity = rms(u);

  % compute standard deviation / RMSD: 
MASvss = single(std(u));
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
temp = single(xcorr(u-MASvsm, rCmax, 'coeff')); 
MASC = single(temp(rCmax + (1:rCmax))); 
fprintf(' done in %.1f seconds.\n', round(10*toc)/10); 

tic

% Find the distance where MASC reaches the correlation value of 1/e.
% We call this distance the oneOverEScale
% We wait till a value less than 1/e has happened twice, just in case
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

matfile = fullfile(fileparts(path), outputname);

%structfile = fullfile(pathname, 'struct.fig');
%histfile = fullfile(pathname, 'hist.fig');
%corelfile = fullfile(pathname, 'corel.fig');
%logcorf = fullfile(pathname, 'lcorel.fig');
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

fprintf('  done in %.1f seconds.\n', round(10*toc)/10);
%close all; 

