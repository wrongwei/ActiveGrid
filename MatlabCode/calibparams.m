% Finds Calibration Parameters for yaw calibration
% must exist already: a1, n1, a2, n2, King's Law coefficients
% found by running loadvelocityff
% Inputs: 
% dirang: angle vector from which calibration was performed (in degrees)
% Volt1, Volt2: wire voltages corresponding to each angle        
% Param (integer) = 1,2,3,4 corresponding to each of the four equations found in Bruun
% 1990-- "Calibration and Analysis of X-hot wire probe signals" equations
% (3(a-d))
% Horace Zhang + Jessie Liu Summer 2014 


function [p1,p2] = calibparams(dirang, Volt1, Volt2, a1, n1, a2, n2, param)
%operational angle (aka alphabar) in radians
opangle = pi/4; 
%convert dirang for each wire to radians and add the operational angle
%these are alpha1,2 in Bruun 
dirangr1 = opangle + dirang.*pi/180 ; 
dirangr2 = opangle - dirang.*pi/180 ;
%find index of 0 degrees, assuming calibration goes from -theta to +theta
ind = round(length(dirang)/2);

%use modified voltage ratio

Volt1s = ((Volt1.^2 - a1)/(Volt1(ind)^2 - a1)).^(1/n1);
Volt2s = ((Volt2.^2 - a2)/(Volt2(ind)^2 - a2)).^(1/n2);

if param == 1
    % Equation 3a in table 1
    y1 = Volt1s.^2-1;
    y2 = Volt2s.^2-1;

    x1 = Volt1s.^2*(sin(opangle))^2 - (sin(dirangr1)).^2;
    x2 = Volt2s.^2*(sin(opangle))^2 - (sin(dirangr2)).^2;
    
    hold on; 
    scatter(x1, y1, 'r', '.');
    scatter(x2, y2, 'b', '.'); 
    legend('wire 1' , 'wire 2');
  
    % find k^2 values
    s1 = polyfit(x1,y1,1);
    s2 = polyfit(x2,y2,1);
    konesq = 1-s1(1);
    ktwosq = 1-s2(1);
    p1 = sqrt(konesq);
    p2 = sqrt(ktwosq);
    
 end;
 
 if param == 2
     %Equation 3b in table 1
     y1 = log(Volt1s);
     y2 = log(Volt2s);
     x1 = log(cos(dirangr1)*1/cos(opangle));
     x2 = log(cos(dirangr2)*1/cos(opangle));
     
     hold on; 
     scatter(x1, y1, 'r', '.');
     scatter(x2, y2, 'b', '.'); 
     legend('wire 1' , 'wire 2');
     %find m
     
     s1 = polyfit(x1,y1,1);
     s2 = polyfit(x2,y2,1);
     p1 = s1(1);
     p2 = s2(1);
 end; 
 
 if param == 3
     %equation 3c in table1
     y1 = sqrt(Volt1s) - 1;
     y2 = sqrt(Volt2s) - 1;
     x1 = sqrt(Volt1s).*(1-sqrt(cos(opangle)))-(1-sqrt(cos(dirangr1)));
     x2 = sqrt(Volt2s).*(1-sqrt(cos(opangle)))-(1-sqrt(cos(dirangr2)));
     
     hold on; 
     scatter(x1, y1, 'r', '.');
     scatter(x2, y2, 'b', '.');
     legend('wire 1' , 'wire 2');
     
     %find b
     s1 = polyfit(x1,y1,1);
     s2 = polyfit(x2,y2,1);
     p1 = s1(1);
     p2 = s2(1);
     
 end;
 
 if param == 4
     %equation 3d in table1
     y1 = cos(dirang.*pi/180) - Volt1s;
     y2 = cos(dirang.*pi/180) - Volt2s;
     x1 = sin(dirang.*pi/180);
     x2 = sin(dirang.*pi/180);
     
     hold on; 
     scatter(x1, y1, 'r', '.');
     scatter(x2, y2, 'b', '.');
     legend('wire 1' , 'wire 2');
     
     %find (theta_e)
     s1 = polyfit(x1,y1,1);
     s2 = polyfit(x2,y2,1);
     %find theta_e 
     p1 = atan(s1(1));
     p2 = atan(s2(1));
     
 end;
 
  
end


