% Script to plot spectrum of velocity data stored in a .mat workspace file
% Written by Nathan Wei and Kevin Griffin, 4 December 2015
% Dependencies: none

path = fileparts('/n/homeserver2/user3a/nwei/Documents/Turbulence2015/data08_24_15/');
addpath(path);

% specify file to load data from
w = load('statscorr_th5.2th50_0824.mat','u','MASvsm','deltaT');

%[MASp, MASf] = hwspectrum(w.u-w.MASvsm, w.deltaT);

samples = length(w.u); 
interval = 5; % seconds of each segment to average

num = samples*w.deltaT/interval;
f = (0:num-1)/(samples-1); 
f = f/w.deltaT; 
spectrum = zeros(num,1);
for i = 1 : num
    tempspectrum = fft(w.u(((i-1)*num+1):(i*num))-w.MASvsm); 
    tempspectrum = real(tempspectrum).^2 + imag(tempspectrum).^2; 
    tempspectrum = tempspectrum * w.deltaT^2; 
    spectrum = spectrum + tempspectrum;
end

% curve fits for sections of spectrum graph
p1 = polyfit(log10(f(2:5)).', log10(spectrum(2:5)), 1);
disp(p1);
p2 = polyfit(log10(f(100:500)).', log10(spectrum(100:500)), 1);
disp(p2);

% plot histogram of velocity fluctuation data
%figure(1);
%histogram(w.u-w.MASvsm);

% plot spectrum
figure(1);
semilogx(f(5:(length(f)-50)), spectrum(5:(length(f)-50)),'.'); % cut off FFT noise
xlabel('Frequency'); % what are the units? Hz? kHz? GigaGriffins?
ylabel('Power');
title('Power Spectrum of Velocity Fluctuations');
