% plots the normalized (by integral length scale) correlation function from data extracted from
% makeallstats.m and saves them in the folder
% Dependencies: hwils.m
% Horace Zhang + Jessie Liu Summer 2014
% Updated by Kevin Griffin and Nathan Wei

% Where is the path?-----------------------------------------------------
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/');
path = fileparts('/Users/kevin/Documents/Data/data08_07_15/');

addpath(path);

% load all the workspaces you want to graph. Put each one in a varaible,
% and then put all of those variables into the array below named
% workspaceArray
% -----------------------------------------------------------------------
close all;
fprintf('Loading workspaces... ');
tic;
%workspace1 = load('statscorr_lt1.3lt1_h0.8_Re_0807.mat');

workspace1 = load('statscorr_lt1.3lt1_h0_Re_0807.mat');
workspace2 = load('statscorr_lt1.3lt1_h0.05_Re_0807.mat');
workspace3 = load('statscorr_lt1.3lt1_h0.1_Re_0807.mat');
workspace4 = load('statscorr_lt1.3lt1_h0.2_Re_0807.mat');
workspace5 = load('statscorr_lt1.3lt1_h0.4_Re_0807.mat');
workspace6 = load('statscorr_lt1.3lt1_h0.8_Re_0807.mat');

workspaceArray = [workspace1, workspace2,workspace3,workspace4,workspace5,workspace6];
workspaceNames = {'rms 30 Height 0','rms 30 Height 0.05','rms 30 Height 0.1',...
    'rms 25 Height 0.2','rms 19 Height 0.4','rms 18 Height 0.8'};
 %workspaceNames = {'works'};   
  chartTitle = 'Correlation Functions for Top Hat Long Tail, SpatialSigma=1.3, TemporalSigma=.1sec, Reynolds Number Constant';

%{
workspace1 = load('statscorr_lt1.3lt1_h0_rms40_0806.mat');
workspace2 = load('statscorr_lt1.3lt1_h0.05_rms40_0806.mat');
workspace3 = load('statscorr_lt1.3lt1_h0.1_rms40_0806.mat');
workspace4 = load('statscorr_lt1.3lt1_h0.2_rms40_0806.mat');
workspace5 = load('statscorr_lt1.3lt1_h0.4_rms40_0806.mat');
workspace6 = load('statscorr_lt1.3lt1_h0.8_rms40_0806.mat');
workspaceArray = [workspace1,workspace2,workspace3,workspace4,workspace5,workspace6];
workspaceNames = {'Height: 0','Height: 0.05',...
    'Height: 0.1','Height: 0.2','Height: 0.4','Height: 0.8'};
chartTitle = 'Correlation Functions for Top Hat Long Tail, SpatialSigma=1.3, TemporalSigma=.1sec, RMS=40deg';
%}

% Example with one workspace
workspace1 = load('statscorr_test_0807.mat');
workspaceArray = [workspace1];
workspaceNames = {'test'};
chartTitle = 'Correlation Function';

%MODIFY THIS, WHAT YOU WANT TO NAME THE
%FIGURES? ---------------------------------------------------------------
figurename = 'lt1.3lt1_rms40_corrfs_lkjhffsdgh.fig';

% This change in involved prefixed the loaded workspace variables with workspaceArray(j).
% These variables include MASC MASvss sepval
% -----------------------------------------------------------------------
if (length(workspaceNames) ~= length(workspaceArray)) % validation
    error('Array size mismatch between workspaceNames and workspaceArray');
end
fprintf('Workspaces loaded in %.1f seconds.\n', round(10*toc)/10);
for j = 1 : length(workspaceArray)
    fprintf('\n-----------------------------------------------------------------------');
    fprintf('------\n Processing test %s \n',workspaceNames{j});
    %sepval = [1:12e6]/(20000)*mean(u); 
    L = hwils(workspaceArray(j).MASC,workspaceArray(j).sepval,2); %this is the integral length scale

    % Cut off MASC (the correlation function) at nth zero crossing. Discards
    % much of unwanted data
    count = 0; 
    n = 2;
    for i=1:length(workspaceArray(j).MASC)
        if(workspaceArray(j).MASC(i) < 0 ) 
            count = count + 1;
        end;
        if(count >= n)
            cutoff = i; 
            break;
        end;

    end;
    MASCc = workspaceArray(j).MASC(1: cutoff); % the cut off structure function
    sepvalc = workspaceArray(j).sepval(1:cutoff); %the cut off sepval 

    nu = 15.11e-6; % kinematic viscosity of air
    % calculate energy dissipation
    epsilon = 0.5 * (workspaceArray(j).MASvss^3) / L; % 0.5 is constant prefactor
    % calculate Kolmogorov length scale
    eta = (nu^0.75) * (epsilon^(-0.25));
    % calculate maximum frequency
    freq = workspaceArray(j).MASvss/eta;
    
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
    %semilogy(sepvalc/L,log(MASCc));
    plot(sepvalc/L,MASCc);
    %plot(sepvalc/taylorL,MASCc);
    %semilogx(sepvalc/eta,MASCc);
    hold on;
    samplePoints = 5:30; % make a vector of the first 26 points excluding noise
    samplePoints2 = sepvalc(samplePoints)/L;
    curvePoints = 1:1000; % vector used for plotting the quadratic fit
    curvePoints = sepvalc(curvePoints)/L;
    corrVals = MASCc(samplePoints); % correlation value of the 26 sample points
    %plot(samplePoints2,corrVals,'-ok'); % plot these 26 points
    hold on;
    p = polyfit(samplePoints2,corrVals',2); %fit a second order polynomial to these 26 points. Note polyfit wanted both vectors to be row vectors, I transpose corrVals
    y1 = polyval(p,curvePoints);
    %plot(curvePoints,y1,'r'); % plot the curve fit
    % print standard deviation = RMS velocity
    fprintf('RMS velocity (m/s) = %f \n', workspaceArray(j).MASvss);
    fprintf('Integral Length Scale = %f\n', L);
    % the x-intercept of polyfit p is the taylor length scale
    taylorL = max(roots(p))*L;
    fprintf('Taylor length scale = %f\n', taylorL);
    %fprintf('Comment: I multipled by L because of the scaling of the graph\n');
    nu = 15.11e-6;
    turbRe = workspaceArray(j).MASvss*taylorL/nu;
    fprintf('Turbulent Reynolds Number = %f\n', turbRe);

    hax = gca; 
    ylabel('correlation   ');
    xlabel('distance (m/L)  ');
    xlim([0 4]);
    ylim([0 1]);
    %title('Correlation Function loglog');
    title(chartTitle);

    %ncorf = fullfile(pathname, 'ncorel_trd1.5.fig');
    logcorf = fullfile(path, figurename);

    %saveas(H1, ncorf);
    saveas(H2, logcorf);
    
    % print results
    fprintf('Energy dissipation rate (W) = %f \n', epsilon);
    fprintf('Kolmogorov length scale (m) = %f \n', eta);
    fprintf('Maximum fluctuation frequency (Hz) = %f \n', freq);
    hold on;
end
    
legend(workspaceNames);
rmpath(path);

% play sound to alert user to end of data processing
t = 0:(1/8000):0.25;
y1 = sin(2*pi*440*t);
y2 = sin(2*pi*880*t);
y = [y1 y2 y1 y1];
sound(y, 8000);
