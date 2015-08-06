
% apply temperature correction to a series of hot wire voltages
%applytempcorrection.m
% 
% Ecorr = applytempcorrection(E, Tm, Tc, Rw, Rc, alpha)
% 
% E - raw wire voltage
% Tm - measured temperature at runtime
% Tc - measured temperature at calibration
% Rw - resistance of wire (hot)
% Rc - resistance of wire (cold)
% alpha - temperature coefficient of resistance (on hot-wire case)
% 
% Greg Bewley  2010
% Revised by Nathan Wei with formulas from Willem van de Water, summer 2015

function Ecorr = applytempcorrection(E, Tm, Tc, Rw, Rc, alpha)
Rcm = Rc * (1 + alpha * (Tm - Tc)); % linear approx. of resistance vs. temp
Ecorr = E * sqrt((Rw - Rc) / (Rw - Rcm)); % multiply voltage by correction
