
% apply temperature correction to a series of hot wire voltages
%applytempcorrection.m
% 
% Ecorr = applytempcorrection(E, T, Twire, T0); 
% 
% E - raw wire voltage
% T - ambient temperature
% Twire - effective wire temperature
%         (the temperature which makes the formula work, 
%          by some calibration)
% T0 - the temperature to correct to.  
% 
% Ecorr = E .* ( (Twire - T0) ./ (Twire - T) ).^(0.5); 
% 
% Greg Bewley  2010

function Ecorr = applytempcorrection(E, T, Twire, T0)

Ecorr = E .* ( (Twire - T0) ./ (Twire - T) ).^(0.5); 
