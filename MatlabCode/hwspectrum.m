
% make power spectrum of hot wire data.  
%hwspectrum.m
% 
% [p freq] = hwspectrum(signal, deltat); 
% 
% 
% Greg Bewley 2009
% 2011:  output is now normalized to have integral 1, 
%  assuming negligible contribution from frequencies above the highest, 
%  rejection of the zero frequency component, 
%  and box-car extrapolation to zero frequency.  

function [p freq] = hwspectrum(signal, deltat)

if nargin < 2, deltat = 1; end

  % way two: 

N = length(signal);   % number of points
T = (N-1)*deltat;     % time interval
freq = (0:N/2-1)/T; 
p = abs(fft(signal))/(N/2); 
p = p(1:N/2).^2; 

specInt = sum(p(2:end).*diff(freq(:))); 
specInt = specInt + freq(2)*p(2); 

  % normalize the integral of the spectrum to 1: 
p = p/specInt; 



% ----
%function [spectrum f p freq] = hwspectrum(signal, deltat)

  % way one: 
%samples = length(signal); 

%f = (0:samples-1)/(samples-1); 
%f = f/deltat; 

%spectrum = fft(signal); 
%spectrum = real(spectrum).^2 + imag(spectrum).^2; 
%spectrum = spectrum * deltat^2; 
