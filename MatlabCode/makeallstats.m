
% script to return standard statistics of a hot wire signal
%makeallstats.m
% 
% must exist: 
% u      - the velocity
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


if (~exist('computestructure')), computestructure = 1; end
if (~exist('highorder')), highorder = 1; end
%fprintf('  working on %d samples.  \n');

%get the directory of your input files:
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/726G0.54/');
pathname = fileparts('/Users/nathan/Documents/Data/data08_03_15/'); % location of calib file
datafolder = fileparts('/Users/nathan/Documents/Data/data08_03_15/g2.6g0.25_10ft/'); % location of data
addpath(pathname);
addpath(datafolder);

%extract velocity
%{
u1 = loadvelocityff('g3g3_0731_00_1.dat', 'calib7_31.m', 1, 1);
u2 = loadvelocityff('g3g3_0731_00_2.dat', 'calib7_31.m', 1, 1);
u3 = loadvelocityff('g3g3_0731_00_3.dat', 'calib7_31.m', 1, 1);
u4 = loadvelocityff('g3g3_0731_00_4.dat', 'calib7_31.m', 1, 1);
%stitch together the file 
u = [u1;u2;u3;u4];
%}

u1 = loadvelocityff('xpos100_ypos100_evts0-2999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u2 = loadvelocityff('xpos100_ypos100_evts3000000-5999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u3 = loadvelocityff('xpos100_ypos100_evts6000000-8999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u4 = loadvelocityff('xpos100_ypos100_evts9000000-11999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u5 = loadvelocityff('xpos100_ypos100_evts12000000-14999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u6 = loadvelocityff('xpos100_ypos100_evts15000000-17999999SN_Ch4.dat', 'calib8_03.m', 1, 1);
u = [u1;u2;u3;u4;u5;u6];

fprintf('velocity extracted \n');

%this is 1/ the sampling frequency
deltaT = 1/20000;

%histX = 35;

% number of correlation function separations starting from one in samples: 
rCmax = 18e6; 

sepval = [1:rCmax]/(20000)*mean(u);
 
  % structure function separations in samples: 
%rs = makelogtime(1, 10000, 50); 
%rs = unique(round(rs)); 

tic
  % compute mean: 
MASvsm = mean(u); 

  % compute rms velocity:
rmsvelocity = rms(u);

  % compute standard deviation: 
MASvss = std(u);
fprintf('  done in %.1f seconds.  spectrum...  \n', round(10*toc)/10); 

  % compute spectrum: 
tic
%[MASp MASf] = hwspectrum(u-MASvsm, deltaT); 
%fprintf('  done in %.1f seconds.  correlation functions...  \n', round(10*toc)/10); 

  % compute correlation functions (as a check, since already computed
  % through S2):
temp = xcorr(u-MASvsm, rCmax, 'coeff'); 
MASC = single(temp(rCmax + (1:rCmax))); 
fprintf('  done in %.1f seconds.  structure functions...  \n', round(10*toc)/10); 

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
fprintf('  done in %.1f seconds.  histogram...  \n', round(10*toc)/10); 

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
fprintf('  done in %.1f seconds.  Saving data...  \n', round(10*toc)/10); 
tic;
matfile = fullfile(pathname, 'statscorr_g2.6g0.25_0803.mat');
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
rmpath(pathname);
rmpath(datafolder); 

fprintf('  done in %.1f seconds.\n', round(10*toc)/10);
%close all; 




