function [ turbRe, taylorL ] = calculatereynoldsnumber(MASCc, sepvalc, MASvss, makePlot)
% take a correlation function named MASC and distance scale named sepval as a arguments and
% calcualte the taylor length scale by fitting a parabola to the first 30
% points of the correlation function. Then use the talyor length scale and
% rms velocity (named MASvss to calcualte reynolds number. Make a plot of
% the parabolic fitting if makePlot is true.

L = hwils(MASCc,sepvalc,2); %this is the integral length scale

%originally it was 5:30 but I prefer 8:50
samplePoints = 8:50; % make a vector of the first 26 points excluding noise
samplePoints2 = sepvalc(samplePoints)/L;
curvePoints = 1:1000; % vector used for plotting the quadratic fit
curvePoints = sepvalc(curvePoints)/L;
corrVals = MASCc(samplePoints); % MASC value of the 26 sample points

p = polyfit(samplePoints2,corrVals',2); %fit a second order polynomial to these 26 points. Note polyfit wanted both vectors to be row vectors, I transpose corrVals
y1 = polyval(p,curvePoints);

% the x-intercept of polyfit p is the taylor length scale
taylorL = max(roots(p)*L);
% multiply by L because of the scaling of the graph
nu = 15.11e-6;
turbRe = MASvss*taylorL/nu;

fprintf('Integral Length Scale = %f\n', L);
fprintf('Taylor length scale = %f\n', taylorL);
fprintf('RMS velocity (m/s) = %f \n', MASvss); % standard deviation = RMS velocity
fprintf('Turbulent Reynolds Number = %f\n', turbRe);

% if L is an argument, then graph the parabolic fitting used to find taylorL
if makePlot
    figure();
    plot(sepvalc/L,MASCc);
    xlabel('distance (m/L)');
    xlim([0 4]);
    ylim([0 1]);
    hold on;
    plot(samplePoints2,corrVals,'-ok'); % plot these 26 points
    plot(curvePoints,y1,'r'); % plot the curve fit
end

end

