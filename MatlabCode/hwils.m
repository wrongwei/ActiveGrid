
% to compute the integral length scale from a correlation function.  
%hwils.m
% function L = hwils(Cr, r, itype); 
% 
% given correlation function C, and its corresponding scale values, 
% return integral length scale L.  
% see code for extrapolation details and parameters.  
% r can be given to specify actual separation values, otherwise 
%  0:length(Cr)-1 will be used, which is usually proportional.  
% itype can be used to specify how the tail end of the correlation 
%  function is contructed: 
%  1: linear fit
%  2: exponential fit (default)
%  3: power law fit
%
%Greg Bewley  7/23/02 revised 2010

function [L, Atail] = hwils(Cr, r, itype)

epsilon = 0.1;  % we will consider the correlation to have reached
                 %  zero when it dips below this value.
%fl = 1/2;	% we will trust data below this fraction of the 
            %  longest length scale...

if nargin < 3
itype = 2;	% 1 for linear fit to construct tail
            % 2 for exponential fit to construct tail
            % 3 for a power law fit.
end

  % for linear fit: 
np = 8;	 % number of points to consider for extrapolation.  
  % for semilogy/powerlaw fit: 
ss = 1/4; % smallest scale we will consider in fit as a fraction of cutoff position.  

if nargin < 2
	r = 0:(length(Cr)-1); 
else
	if isempty(r)
		r = 0:(length(Cr)-1); 
	else	
		[r, I] = sort(r); 
		r = [0 r(:)']; 
		Cr = [1 Cr(I)']; 
	end
end

if nargout == 0
	plotme = 1; 
else
	plotme = 0; 
end

np = np-1; 

  % find fraction fl of the largest length scale, and matrix location:
%lr = fl*max(r);
%lc = find(r>lr);
%lc = min(lc);
lc = length(r); 

  % find where the function dips below epsilon: 
bec = find(Cr<epsilon); 
  % choose limit of integration: 
if isempty(bec),      bec = lc; 
%elseif min(bec)>lc,  bec = lc; 
else,                 bec = min(bec); 
end
ss = round(bec*ss); 

%bec, ss

if plotme
	figure, plot(r,Cr, 'x--')
	hold on
end

if itype == 1
	  % compute a linear fit to bring the correlation to zero:
	if bec <= np, error('not enough data to extrapolate.'), end
	pe = polyfit(r(bec-np:bec),Cr(bec-np:bec),1);
	r0 = -pe(2)/pe(1);
	  % find tail area:
	Atail = (Cr(bec)/2)*(r0-r(bec));
	if plotme
		plot([r(bec), r0], [Cr(bec), 0], 'r+-')
	end
elseif itype == 2
	  % compute semilogy fit to bring correlation to zero: 
	if bec <=ss, error('not enough data to extrapolate.'), end
	pe = polyfit(r(ss:bec), log(Cr(ss:bec)), 1); 
	  % find tail area:	
	Atail = -exp(pe(2))*exp(r(bec)*pe(1))/pe(1); 
	  % extend r for plot: 
	r0 = log(0.001/exp(pe(2)))/pe(1); 
	rp = r(ss):(r(2)-r(1)):r0; 
	if plotme
		plot(rp, exp(pe(2))*exp(pe(1)*rp), 'r--')
		plot(r(bec), exp(pe(2))*exp(pe(1)*r(bec)), 'r+')
	end
elseif itype == 3
	  % compute power law fit: 
	if bec <= ss, error('not enough data to extrapolate.'), end
	pe = findpowerlaw(r(ss:bec), Cr(ss:bec)); 
	  % find tail area: 
	Atail = (pe(2)/(pe(1)+1)) * r(bec)^(pe(1)+1); 
	  % extend r for plot: 
	r0 = (0.05/pe(2))^(1/pe(1)); 
	rp = r(ss):(r(2)-r(1))*100:r0; 
	if plotme
		plot(rp, pe(2)*rp.^pe(1),'r--')
		plot(r(bec), pe(2)*r(bec)^pe(1),'r+')
	end
end

  % compute integral of bulk
Cint = 0; 
for i = 1:bec-1
	dr = r(i+1)-r(i); 
	aCr = (Cr(i)+Cr(i+1))/2; 
	Cint = Cint + aCr*dr; 
end

  % length scale is total area: 
L = Cint + Atail; 
