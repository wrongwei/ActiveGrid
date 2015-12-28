% Script to plot spectrum of velocity data stored in a .mat workspace file
% Written by Nathan Wei and Kevin Griffin, 4 December 2015
% Dependencies: none

path = fileparts('/n/homeserver2/user3a/nwei/Documents/Turbulence2015/data08_24_15/');
addpath(path);

% specify file to load data from
w = load('statscorr_th5.2th50_0824.mat','u','MASvsm','deltaT');

%[MASp, MASf] = hwspectrum(w.u-w.MASvsm, w.deltaT);

samples = length(w.u); 

num = samples*w.deltaT/10;
f = (0:num-1)/(samples-1); 
f = f/w.deltaT; 
spectrum = zeros(num,1);
for i = 1 : num
    tempspectrum = fft(w.u(((i-1)*num+1):(i*num))-w.MASvsm); 
    tempspectrum = real(tempspectrum).^2 + imag(tempspectrum).^2; 
    tempspectrum = tempspectrum * w.deltaT^2; 
    spectrum = spectrum + tempspectrum;
end

% plot histogram of velocity fluctuation data
%figure(1);
%histogram(w.u-w.MASvsm);

% plot spectrum
figure(1);
semilogx(f, spectrum);
xlim([0 max(f)*0.95]); % cut off noise from FFT
xlabel('Frequency'); % what are the units? Hz? kHz? GigaGriffins?
ylabel('Power');
title('Power Spectrum of Velocity Fluctuations');
