
% calibration file processing program
%processcalibrationff.m
% 
% SINGLE CHANNEL VERSION OF PROCESSCALIBRATIONNF: 
% 
% [a b n R figs] = processcalibrationff(calibrationfilename, T0, freeexponent, ...  
%                                       probenumber, sensornumber, figs); 
% 
% T0 (optional)   - temperature to correct each data point to.  
%                   default no temp correction.  
% freeexponent (optional) - value for King's law exponent, or true if left free to fit.  
%                           default 0.5.  
% figs (optional) - 0 to make one new figure for each probe, 
%                   or a figure handle to add to for each probe.  
%                   default no figs.  
% 
% Gregory Bewley  2009

function [a b n R figs] = processcalibrationff(calibrationfilename, T0, freeexponent, ...  
                                               probenumber, sensornumber, figs)

if (nargin < 6), figs = 0; end
%if (nargin < 3), freeexponent = 0.5; 
%elseif isempty(freeexponent), freeexponent = 0.5; end
if (nargin < 3), freeexponent = true; 
elseif isempty(freeexponent), freeexponent = 0.5; end
if (nargin < 2), T0 = []; end

  % ---- load probeinfo, calibdata, ...  
[pathstr, name] = fileparts(calibrationfilename); 
odiry = pwd; 
if (~isempty(pathstr)), cd(pathstr); end
run(name); 
cd(odiry); 

  % ---- organize the data: 
  % maximum channel number: 
maxchannel = 0; 
for i = 1:probeinfo.numprobes
	for j = 1:probeinfo.probe(i).numsensors
		if (probeinfo.probe(i).sensor(j).channel > maxchannel)
			maxchannel = probeinfo.probe(i).sensor(j).channel; 
		end
	end
end

if (size(calibdata.v, 1) == 1)
	calibdata.v = repmat(calibdata.v, [maxchannel 1]); 
end
if (size(calibdata.T, 1) == 1)
	calibdata.T = repmat(calibdata.T, [maxchannel 1]); 
end

  % ---- sort data by increasing velocity: 
channel = probeinfo.probe(probenumber).sensor(sensornumber).channel; 
[calibdata.v(channel, :), I] = sort(calibdata.v(channel, :)); 
if ~isempty(calibdata.fanspeed),   calibdata.fanspeed = calibdata.fanspeed(I);   end
if ~isempty(calibdata.dp),   calibdata.dp = calibdata.dp(I);   end
if ~isempty(calibdata.dens) && length(calibdata.dens) > 1
	calibdata.dens = calibdata.dens(I); 
end
if (length(calibdata.T(channel, :)) > 1)
	calibdata.T(channel, :) = calibdata.T(channel, I); 
end
calibdata.E(channel, :) = calibdata.E(channel, I); 

  % the relevant data: 
thisE = calibdata.E(channel, :); 

  % ---- apply temperature correction, if requested: 
if ~isempty(T0)
    calibdata.T(channel, :); 
	thisE = applytempcorrection(thisE, calibdata.T(channel, :), calibdata.Twire(channel), T0); 
end

  % ---- find King's law coefficients for the channel: 
a = []; b = []; n = []; R = []; 
[a b n Rtemp] = findKingslaw(calibdata.v(channel, :), thisE, freeexponent); 
%R = mean(std(Rtemp-thisE)./thisE); %this line of code causes errors

  % ---- plot the data, if requested: 
if (nargin > 5)
	  % make figures
	  % one figure for each probe
	if (figs == 0)
			figs = figure; 
            fig2 = figure; 
	end
	
	figure(figs); 
	title(probeinfo.probe(probenumber).probename); 
	hold on
    
	figure(fig2); 
	title(probeinfo.probe(probenumber).probename); 
	hold on
    
    E = calibdata.E(channel, :); 
      % ---- apply temperature correction, if requested: 
    if ~isempty(T0)
        E = applytempcorrection(E, calibdata.T(channel, :), calibdata.Twire(channel), T0); 
    end

    u = calibdata.v(channel,:); 
    meandens = mean(calibdata.dens); 

%			plot(u, E, 'x', ...  
%			plot(u, ((E.^2 - a)/b).^2, 'x', ...  
    figure(fig2)
    plot(u, E.^2, 'bx', 'MarkerSize', 16)
    plot(u, a + b*u.^n, 'k-')
    
    figure(figs)
    plot(u, (E.^2 - a).^(1/n)/meandens, 'x', 'MarkerSize', 16); 

%			plotKingslaw([u(1)*0.9 u(end)*1.1], a(channel), b(channel), n(channel), figs(i)); 
    us = u(1)*0.9:(u(end)*1.1-u(1)*0.9)/100:u(end)*1.1; 
    plot(us, us*(b).^(1/n)/meandens, 'k-'); 

    h = xlabel('u [m/s]'); 
    setfontsize(h, 18); 
    h = ylabel('(E^2 - a)^{1/ n} / \rho ,  u b^{1/ n} / \rho'); 
    setfontsize(h, 18); 
    setfontsize(gca, 16); 
end
