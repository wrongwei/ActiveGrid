% plots the normalized (by integral length scale) correlation function from data extracted from
% makeallstats.m and saves them in the folder
% Dependencies: hwils.m
% Horace Zhang + Jessie Liu Summer 2014

% Where is the path?-----------------------------------------------------
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/');
path = fileparts('/Users/nathan/Documents/Data/data08_05_15/');
addpath(path); 

% which workspace stats do you want to load??????????
% -----------------------------------------------------------------------
load statscorr_g1g1_0805.mat
% -----------------------------------------------------------------------

%sepval = [1:12e6]/(20000)*mean(u); 
L = hwils(MASC,sepval,2) %this is the integral length scale

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
plot(sepvalc/L,MASCc);
hax = gca; 
ylabel('correlation   ');
xlabel('distance (m/L)  ');
xlim([0 4]);
%title('Correlation Function loglog');
title('Correlation Function');

%MODIFY THIS, WHAT YOU WANT TO NAME THE
%FIGURES? -------------------------------------------------------------------------------
%ncorf = fullfile(pathname, 'ncorel_trd1.5.fig');
logcorf = fullfile(path, 'test.fig');

%saveas(H1, ncorf);
saveas(H2, logcorf);


rmpath(path);
 
