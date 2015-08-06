% plots the normalized (by integral length scale) correlation function from data extracted from
% makeallstats.m and saves them in the folder
% Dependencies: hwils.m
% Horace Zhang + Jessie Liu Summer 2014
% Updated by Kevin Griffin and Nathan Wei

% Where is the path?-----------------------------------------------------
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/');
path = fileparts('/Users/kevin/Documents/Data/data07_31_15/');
addpath(path);
path = fileparts('/Users/kevin/Documents/Data/data08_03_15/');
addpath(path);
% load all the workspaces you want to graph. Put each one in a varaible,
% and then put all of those variables into the array below named
% workspaceArray
% -----------------------------------------------------------------------
close all;
fprintf('Loading workspaces... ');
tic;
% Example with three workspaces
workspace1 = load('statscorr_g2.6g1_0803.mat');
workspace2 = load('statscorr_g2.6inv2_0731.mat');
workspace3 = load('statscorr_g2.6th1_0731.mat');
workspace4 = load('statscorr_g2.6tr2_0731.mat');
workspace5 = load('statscorr_randg1_0731.mat');
workspace6 = load('statscorr_us2.6g1_0731.mat');
workspaceArray = [workspace1,workspace2,workspace3,workspace4,workspace5,workspace6];
workspaceNames = {'g2.6g1','g2.6inv2','g2.6th1',...
    'g2.6tr2','randg1','us2.6g1'};
chartTitle = 'Correlation Functions 0731 and 0803';
% Example with one workspace
%workspace1 = load('statscorr_g1g1_0805.mat');
%workspaceArray = [workspace1];

%MODIFY THIS, WHAT YOU WANT TO NAME THE
%FIGURES? ---------------------------------------------------------------
figurename = 'g0.5_corrfs.fig';

% This change in involved prefixed the loaded workspace variables with workspaceArray(j).
% These variables include MASC MASvss sepval
% -----------------------------------------------------------------------
if (length(workspaceNames) ~= length(workspaceArray)) % validation
    error('Array size mismatch between workspaceNames and workspaceArray');
end
fprintf('Workspaces loaded in %.1f seconds.\n', round(10*toc)/10);
for j = 1 : length(workspaceArray)
    fprintf('\n Processing test %s \n',workspaceNames{j});
    %sepval = [1:12e6]/(20000)*mean(u); 
    L = hwils(workspaceArray(j).MASC,workspaceArray(j).sepval,2); %this is the integral length scale

    nu = 1.46e-5; % kinematic viscosity of air
    % print standard deviation = RMS velocity
    fprintf('RMS velocity (m/s) = %.4f \n', workspaceArray(j).MASvss);
    % calculate energy dissipation
    epsilon = 0.5 * (workspaceArray(j).MASvss^3) / L; % 0.5 is constant prefactor
    % calculate Kolmogorov length scale
    eta = (nu^0.75) * (epsilon^(-0.25));
    % calculate maximum frequency
    freq = workspaceArray(j).MASvss/eta;
    % print results
    fprintf('Integral length scale (m) = %.8f \n', L);
    fprintf('Energy dissipation rate (W) = %.8f \n', epsilon);
    fprintf('Kolmogorov length scale (m) = %.8f \n', eta);
    fprintf('Maximum fluctuation frequency (Hz) = %.8f \n', freq);
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
    hold on;
    samplePoints = 5:30; % make a vector of the first 26 points excluding noise
    samplePoints2 = sepvalc(samplePoints)/L;
    curvePoints = 1:1000; % vector used for plotting the quadratic fit
    curvePoints = sepvalc(curvePoints)/L;
    corrVals = MASCc(samplePoints); % correlation value of the 26 sample points
    plot(samplePoints2,corrVals,'-ok'); % plot these 26 points
    hold on;
    p = polyfit(samplePoints2,corrVals',2); %fit a second order polynomial to these 26 points. Note polyfit wanted both vectors to be row vectors, I transpose corrVals
    y1 = polyval(p,curvePoints);
    plot(curvePoints,y1,'r'); % plot the curve fit
    fprintf('Integral Length Scale = %f\n', L);
    % the x-intercept of polyfit p is the taylor length scale
    taylorL = max(roots(p))*L;
    fprintf('The Taylor length scale = %f\n', taylorL);
    fprintf('Comment: I multipled by L because of the scaling of the graph\n');
    turbRe = workspaceArray(j).MASvss*taylorL/nu;
    fprintf('mu = %f\nrms = %f\nTurbulent Reynolds Number = %f\n',nu, workspaceArray(j).MASvss, turbRe);

    hax = gca; 
    ylabel('correlation   ');
    xlabel('distance (m/L)  ');
    if (j == length(workspaceArray))
        %legend(workspaceNames);
    end
    xlim([0 4]);
    ylim([0 1]);
    %title('Correlation Function loglog');
    title(chartTitle);

    %ncorf = fullfile(pathname, 'ncorel_trd1.5.fig');
    logcorf = fullfile(path, figurename);

    %saveas(H1, ncorf);
    saveas(H2, logcorf);
end

rmpath(path);
