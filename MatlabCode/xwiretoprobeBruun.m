
% convert X-wire hot wire velocities to probe coordinate velocities.  
%xwiretoprobeBruun.m
% 
% [vs vn] = xwiretoprobeBruun(v1, v2, alpha1, alpha2); 
% 
% from Bruun (1996) Meas. Sci. Technol.  
% see also Bruun et al. (1990) Meas. Sci. Technol.  
% 
% Greg Bewley  2010
% alpha1 and alpha2 are the operational angles of the wires (alphabar)
% in Bruun 1990 (Calibration and analysis of X hot wire probe signals)

function [vs vn] = xwiretoprobeBruun(v1, v2, alpha1, alpha2)

  % the yaw factors for the wires: 
k1 = sqrt(0.075); 
%k1 = 0; 
k2 = sqrt(0.017); 
%k2 = 0; 

  % the wire angles: 
%if nargin < 4, alpha2 = pi/4; end
%if nargin < 3, alpha1 = pi/4; end

  % convert degrees to radians, if necessary: 
if (alpha1 > 2), alpha1 = alpha1*pi/180; end
if (alpha2 > 2), alpha2 = alpha2*pi/180; end

  % coefficients, angular response functions: 
  % from Bruun (1996) Meas. Sci. Technol.  
  % see also Bruun et al. (1990) Meas. Sci. Technol.  
f1 = sqrt( (cos(alpha1))^2 + k1^2 * (sin(alpha1))^2 ); 
f2 = sqrt( (cos(alpha2))^2 + k2^2 * (sin(alpha2))^2 ); 
%g1 = (1-k1^2)*(cos(alpha1))^2 * tan(alpha1) / ((cos(alpha1))^2 + k1^2 * (sin(alpha1))^2); 
%g2 = (1-k2^2)*(cos(alpha2))^2 * tan(alpha2) / ((cos(alpha2))^2 + k2^2 * (sin(alpha2))^2); 
%same as above but in a different form
g1 = (1-k1^2)*(cos(alpha1))^2 * tan(alpha1) / ((cos(alpha1))^2 * (1 - k1^2) + k1^2); 
g2 = (1-k2^2)*(cos(alpha2))^2 * tan(alpha2) / ((cos(alpha2))^2 * (1 - k2^2) + k2^2); 

  % effective cooling velocities: 
v1e = v1 * f1; 
v2e = v2 * f2;  
   %(see Bruun 1990 equation 2) These are exactly the effective cooling
   % velocities

  % the velocities in the probe coordinates: 
%vs = ( g2*v1e/f1 + g1*v2e/f2 ) / (g1 + g2); 
%vn = (    v2e/f2 -    v1e/f1 ) / (g1 + g2); 
 % correction March 5th 2013: cooling velocity is the one measured!  
 % See Bruun equations (10a) and (10b)
vs = ( g2*v1/f1 + g1*v2/f2 ) / (g1 + g2); 
vn = (    v2/f2 -    v1/f1 ) / (g1 + g2); 

