
% load in velocities from data file, check isotropy using user choice method
% 
% [iso] = checkIso(pn, inp)
% 
% from Bruun (1996) Meas. Sci. Technol.  
% see also Bruun et al. (1990) Meas. Sci. Technol.  
% 
% Greg Bewley  2010
% alpha1 and alpha2 are the operational angles of the wires (alphabar)
% in Bruun 1990 (Calibration and analysis of X hot wire probe signals)
% Modified Horace Zhang + Jessie Liu Summer 2014
%
% Dependencies: calibparams.m, loadvelocityff.m 
% Inputs are
% pn (String): the pathname to data file 
% inp (int): which method do you want to use?
% output is std(streamwisevelocity)/std(normalvelocity) should be = 1 in
% the isotropic case
% NOTE: a graph also pops up, don't worry about it unless it is
% very nonlinear
% Keep in mind which channel is which wire (don't reverse them!)
% You will know that you reversed them if the graphs look (very) non-linear
% and/or the graphs for inp = 1,2 or 3 will have negative slope
% must manually enter calibration angles and voltages

function [iso] = checkIso(pn, inp)
 pn = strcat(pn,'/');
 pathname = fileparts(pn);
 addpath(pathname) 
 % load in velocities-------------------------------------------------------
 [u1 a1 b1 n1] = loadvelocityff('xpos100_ypos100_evts0-2999999XCh2_Ch4.dat', 'Calibxwire7_18.m', 1, 2); 
 u2 = loadvelocityff('xpos100_ypos100_evts3000000-5999999XCh2_Ch4.dat', 'Calibxwire7_18.m', 1, 2);
 u3 = loadvelocityff('xpos100_ypos100_evts6000000-8999999XCh2_Ch4.dat', 'Calibxwire7_18.m', 1, 2);
 u4 = loadvelocityff('xpos100_ypos100_evts9000000-11999999XCh2_Ch4.dat', 'Calibxwire7_18.m', 1, 2);
 u = [u1;u2;u3;u4];

 [v1 a2 b2 n2]  = loadvelocityff('xpos100_ypos100_evts0-2999999XCh1_Ch1.dat', 'Calibxwire7_18.m', 1, 1);
 v2 = loadvelocityff('xpos100_ypos100_evts3000000-5999999XCh1_Ch1.dat', 'Calibxwire7_18.m', 1, 1);
 v3 = loadvelocityff('xpos100_ypos100_evts6000000-8999999XCh1_Ch1.dat', 'Calibxwire7_18.m', 1, 1);
 v4 = loadvelocityff('xpos100_ypos100_evts9000000-11999999XCh1_Ch1.dat', 'Calibxwire7_18.m', 1, 1);
 v = [v1;v2;v3;v4];
 rmpath(pathname);
 
 %Calibration angles (manually entered) -------------------------------------------
 dirang = [-15 -10 -5 0 5 10 15];
 
 %Voltages for calibration angles (manually entered)-----------------------------------
 
 dv1 = [1.769 1.760 1.750 1.738 1.724 1.713 1.700];
 dv2 = [1.626 1.638 1.649 1.661 1.672 1.687 1.689];
 % ------------------------------------------------------------------------------------
 %the wire angles

 alpha1 = pi/4;
 alpha2 = pi/4;

if inp == 1
  % Equation 3a from Table 3 Buun 1990
  % coefficients, angular response functions: 
  % from Bruun (1996) Meas. Sci. Technol.  
  % see also Bruun et al. (1990) Meas. Sci. Technol.  
    [k1 k2] = calibparams(dirang,dv1,dv2,a1,n1,a2,n2,1);
    f1 = sqrt( (cos(alpha1))^2 + k1^2 * (sin(alpha1))^2 ); 
    f2 = sqrt( (cos(alpha2))^2 + k2^2 * (sin(alpha2))^2 ); 

    g1 = (1-k1^2)*(cos(alpha1))^2 * tan(alpha1) / ((cos(alpha1))^2 * (1 - k1^2) + k1^2); 
    g2 = (1-k2^2)*(cos(alpha2))^2 * tan(alpha2) / ((cos(alpha2))^2 * (1 - k2^2) + k2^2); 

    % effective cooling velocities: 
    %v1e = u * f1; 
    %v2e = v * f2;  
    % correction March 5th 2013: cooling velocity is the one measured!  
    % See Bruun equations (10a) and (10b)
    
    vs = ( g2*u/f1 + g1*v/f2 ) / (g1 + g2); 
    vn = (    v/f2 -    u/f1 ) / (g1 + g2);
    std(vs)
    std(vn)
    iso = std(vs)/std(vn)
    
    
end;

if inp == 2
    %Equation 3b
    [m1 m2] = calibparams(dirang,dv1,dv2,a1,n1,a2,n2,2);
    f1 = cos(alpha1)^(m1);
    f2 = cos(alpha2)^(m2);
    g1 = m1*tan(alpha1);
    g2 = m2*tan(alpha2);
    vs = ( g2*u/f1 + g1*v/f2 ) / (g1 + g2); 
    vn = (    v/f2 -    u/f1 ) / (g1 + g2); 
    iso = std(vs)/std(vn)
end;

if inp == 3
    %Equation 3c
    [b1 b2] = calibparams(dirang,dv1,dv2,a1,n1,a2,n2,3);
    f1 = sqrt((1-b1*(1- cos(alpha1))));
    f2 = sqrt((1-b2*(1- cos(alpha2))));
    g1 = (b1*(1-b1)*sqrt(cos(alpha1))+b1^2*cos(alpha1))/(1-b1*(1-sqrt(cos(alpha1))))^2*tan(alpha1);
    g2 = (b2*(1-b2)*sqrt(cos(alpha2))+b1^2*cos(alpha2))/(1-b2*(1-sqrt(cos(alpha2))))^2*tan(alpha2);
    vs = ( g2*u/f1 + g1*v/f2 ) / (g1 + g2); 
    vn = (    v/f2 -    u/f1 ) / (g1 + g2); 
    iso = std(vs)/std(vn)
end;

if inp == 4
    %Equation 3d
    [ae1 ae2] = calibparams(dirang,dv1,dv2,a1,n1,a2,n2,4);
    f1 = cos(ae1);
    f2 = cos(ae2);
    g1 = tan(ae1);
    g2 = tan(ae2);
    vs = ( g2*u/f1 + g1*v/f2 ) / (g1 + g2); 
    vn = (    v/f2 -    u/f1 ) / (g1 + g2); 
    iso = std(vs)/std(vn)
end;

end
