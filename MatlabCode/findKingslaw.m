
% to find King's law from a hot wire calibration data series.  
%findKingslaw.m
% 
% [a b n] = findKingslaw(u, E, freeexponent); 
% 
% fit: 
% E^2 = a + b * u ^ n
% 
% freeexponent (optional): 
%  if true:  exponent is free parameter in fit.  
%  if false: exponent fixed to 0.5
%  if 0:     exponent is free parameter in fit.  
%  if n:     exponent is fixed to n
%  defaults to exponent as free parameter in fit.  
% 
% Gregory Bewley  2009

function [a b n R] = findKingslaw(u, E, freeexponent)

  % this complicated front end is to maintain backward 
  %  compatibility: 
if nargin < 3, freeexponent = true; 
elseif isa(freeexponent, 'logical')
	if (freeexponent == false), exponent = 0.5;  end
else
	if (freeexponent == 0), freeexponent = true; 
	else, exponent = freeexponent; freeexponent = false; end
end

options = statset('MaxIter', 10000, ...  
                  'TolX', 1e-4, 'TolTypeX', 'abs', ...  
                  'TolFun', 1e-5, 'TolTypeFun', 'abs'); 

  % try nlinfit: 
if freeexponent
	x0 = [0.1 0.1 0.5]; 
	[x, R] = nlinfit(u, E, @kingsfunctionfree, x0, options); 
	a = x(1); 
	b = x(2); 
	n = x(3); 
else
exponent
	x0 = [0.1 0.1]; 
	[x, R] = nlinfit(u, E, @(x, u) kingsfunctionfixed(x, u, exponent), x0, options); 
	a = x(1); 
	b = x(2); 
	n = exponent; 
end


% --
% gives the fitted function values

function yfit = kingsfunctionfree(x, u)

yfit = x(1) + x(2)*u.^x(3); 
yfit = (yfit).^(1/2); 


% --
% gives the fitted function values

function yfit = kingsfunctionfixed(x, u, exponent)

yfit = x(1) + x(2)*u.^exponent; 
yfit = (yfit).^(1/2); 



%fprintf('%f %f %f\n', x(1), x(2), x(3)); 
%plot(u, x(1) + x(2)*u.^x(3), 'x')
%pause

% ---- old attempts: 
%c = [u, E]; 

  % fminsearch doesn't seem to work: 
%x = fminsearch(@(x) kingsfunctional(x, c), x0); 

% --
% function gives the mean difference of squares

function Delta = kingsfunctional(x, c)

u = c(1:end/2); 
E = c(end/2 + 1 : end); 
Delta = sum( E.^2 - (x(1) + x(2)*u.^x(3)) ); 
fprintf('%f %f %f %f\n', x(1), x(2), x(3), Delta); 
plot(u, x(1) + x(2)*u.^x(3), 'x')
pause
