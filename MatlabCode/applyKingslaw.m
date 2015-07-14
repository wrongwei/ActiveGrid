
% apply King's law to transform hot wire voltages to velocities.  
%applyKingslaw.m
% 
% u = applyKingslaw(E, a, b, n); 
% 
% Gregory Bewley  2009

function u = applyKingslaw(E, a, b, n)

u = ((E.^2 - a) / b).^(1/n); 
