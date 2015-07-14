
% spectrum smoothing program
%smoothspectrumnew.m
% 
% [spectrumfil, freqfil] = smoothspectrumnew(spectrum, freq, numsamples, overlap); 
% 
% numsamples logarithmically spaced bins spaced between left and right 
%   extremes of frequencies freq.  
% 
% optional overlap specifies the amount of overlap between adjacent bins 
%   between 0 and 1.  
% 
% Greg Bewley 2011

function [spectrumfil, freqfil] = smoothspectrumnew(spectrum, freq, numsamples, overlap)

if nargin < 4, overlap = 0;      end
if nargin < 3, numsamples = 200; end

  % frequency edges: 
edges = unique(round(makelogtime(1, length(freq), numsamples+1))); 

speclength = length(spectrum); 

spectrumfil(1) = spectrum(1); 
freqfil(1) = freq(1); 
for i = 2:length(edges)-1
	
	spectrumfil(i) = mean(spectrum(edges(i):edges(i+1))); 
	freqfil(i) = mean(freq(edges(i):edges(i+1))); 
	
end
