% Script to plot spectrum of velocity data stored in a .mat workspace file
% Written by Nathan Wei and Kevin Griffin, 4 December 2015
% Dependencies: hwspectrum.m

path = fileparts('/n/homeserver2/user3a/nwei/Documents/Turbulence2015/data08_24_15/');
addpath(path);

% specify file to load data from
w = load('statscorr_th5.2th50_0824.mat','u','MASvsm','deltaT');

[MASp, MASf] = hwspectrum(w.u-w.MASvsm, w.deltaT);

% plot histogram of velocity fluctuation data
figure(1);
histogram(w.u-w.MASvsm);

% plot spectrum
figure(2);
plot(MASf, MASp);
xlim([0 0.4]);
xlabel('Frequency'); % what are the units? Hz? kHz? GigaGriffins?
ylabel('Power');
title('Power Spectrum of Velocity Fluctuations');
