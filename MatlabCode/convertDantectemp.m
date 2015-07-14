
% convert Dantec temperature in raw volts to temperature in degrees C
%convertDantectemp.m
% 
% T = convertDantectemp(Tvolts); 
% 
% Greg Bewley  2010

function T = convertDantectemp(Tvolts)

  % use numbers taken from the Dantec Streamware software, 
  % with the internal conversion turned off, then on at various temperatures: 
Vdata = [0.5513 0.6106 0.8445 0.7471 0.7078]; 
Tdata = [16.54  18.31  25.34  22.41  21.23 ]; 

p = polyfit(Vdata, Tdata, 1); 

T = Tvolts*p(1) + p(2); 
