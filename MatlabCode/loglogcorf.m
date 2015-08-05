% plots the normalized (by integral length scale) correlation function from data extracted from
% makeallstats.m and saves them in the folder
% Dependencies: hwils.m
% Horace Zhang + Jessie Liu Summer 2014

% Where is the path?-----------------------------------------------------
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/');
path = fileparts('/Users/kevin/Documents/Data/data08_05_15/');
addpath(path); 

% which workspace stats do you want to load??????????
% -----------------------------------------------------------------------
load statscorr_g1g1_0805.mat
% -----------------------------------------------------------------------

%sepval = [1:12e6]/(20000)*mean(u); 
L = hwils(MASC,sepval,2); %this is the integral length scale

% Cut off MASC (the correlation function) at nth zero crossing. Discards
% much of unwanted data
count = 0; 
n = 2;
for i=1:length(MASC)
    if(MASC(i) < 0 ) 
        count = count + 1;
    end;
    if(count >= n)
        cutoff = i; 
        break;
    end;
    
end;
MASCc = MASC(1: cutoff); % the cut off structure function
sepvalc = sepval(1:cutoff); %the cut off sepval 

% Plot correlation function on normal/loglog axes. with x-axis normalized
% by integral lengthscale
% gives a vertical line on where the integral length scale is on the
% correlation function.

hold on;
%{
H1 = figure(1);
grid on;
set(gca, 'fontsize', 15);
plot(sepvalc/L,MASCc, 'LineWidth', 2);
hax = gca; 
line([L,L],get(hax,'YLim'), 'Color' , [0 0 0], 'LineWidth', 1); 
ylabel('correlation   ');
xlabel('distance (m/L)  ');
title('Correlation Function')
%}

H2 = figure(2);
set(gca, 'fontsize', 12);
%loglog(sepvalc/L,MASCc,'*');
plot(sepvalc/L,MASCc,'b');
hold on;
distR = 5:30; % make a vector of the first 26 points excluding noise
linux = 1:1000;
linux2 = sepvalc(linux)/L;
distR2 = sepvalc(distR)/L;
corrVals = MASCc(distR);
corrVals = corrVals';
plot(distR2,corrVals,'-ok'); %plot these 26 points
hold on;
p = polyfit(distR2,corrVals,2); %fit a second order polynomial to these 26 points
y1 = polyval(p,linux2);
plot(linux2,y1,'r');
fprintf('Integral Length Scale = %f\n', L);
% the x-intercept of polyfit p is the taylor length scale
taylorL = max(roots(p))*L;
fprintf('The taylor length scale = %f\n', taylorL);
fprintf('Comment: I multipled by L because of the scaling of the graph\n');
mu = 1.46E-5;
turbRe = MASvss*taylorL/mu;
fprintf('mu = %f\nrms = %f\nTurbulent reynolds Number = %f\n',mu, MASvss, turbRe);

hax = gca; 
ylabel('correlation   ');
xlabel('distance (m/L)  ');
xlim([0 4]);
ylim([0 1]);
%title('Correlation Function loglog');
title('Correlation Function');

%MODIFY THIS, WHAT YOU WANT TO NAME THE
%FIGURES? -------------------------------------------------------------------------------
%ncorf = fullfile(pathname, 'ncorel_trd1.5.fig');
logcorf = fullfile(path, 'test.fig');

%saveas(H1, ncorf);
saveas(H2, logcorf);


rmpath(path);
 
