% plots the normalized (by integral length scale) correlation function from data extracted from
% makeallstats.m and saves them in the folder
% Dependencies: hwils.m
% Horace Zhang + Jessie Liu Summer 2014
% Updated by Kevin Griffin and Nathan Wei

% Where is the path?-----------------------------------------------------
%pathname = fileparts('/Users/Horace/Documents/Germany2014/MATLABCode/MoreCode/DecayData/');

path = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_05_15/');
addpath(path);
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_06_15/');
addpath(path1);
path2 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_07_15/');
addpath(path2);
path3 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_10_15/');
addpath(path3);
path4 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_11_15/');
addpath(path4);
path5 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_12_15/');
addpath(path5);
path6 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_17_15/');
addpath(path6);
path7 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_18_15/');
addpath(path7);
% load all the workspaces you want to graph. Put each one in a varaible,
% and then put all of those variables into the array below named
% workspaceArray
% -----------------------------------------------------------------------

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('*****************************************************************************\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Loading workspaces... ');
tic;

%workspace2 = load('statscorr_lt5.2lt4_0812.mat','MASC','sepval','MASvss','oneOverEScale');
%workspace2 = load('statscorr_lt5.2lt50_0812.mat','MASC','sepval','MASvss','oneOverEScale');
%workspace3 = load('statscorr_th2.6th2_0812.mat','MASC','sepval','MASvss','oneOverEScale');

workspace2 = load('statscorr_th3.9th3_0817.mat','MASC','sepval','MASvss','oneOverEScale');
%{
workspace5 = load('statscorr_us3us2_0812.mat','MASC','sepval','MASvss','oneOverEScale');
workspace6 = load('statscorr_us3us2_oldimp_newus_0818.mat','MASC','sepval','MASvss','oneOverEScale');
workspace7 = load('statscorr_th3th2_oldimp_abs_0818.mat','MASC','sepval','MASvss','oneOverEScale');
workspace8 = load('statscorr_us3us2_oldimp_abs_0818.mat','MASC','sepval','MASvss','oneOverEScale');
%}
workspaceArray = [...
    
    %workspace1...
    %,...
    workspace2...
    
    %,...
    %workspace3...
    %{
    ,...
    workspace4...
    ,...
    workspace5...
    ,...
    workspace6...
    ,...
    workspace7...
    ,...
    workspace8...
    %}
    ];

workspaceNames = {...
    %'lt5.2lt4',...
    'lt5.2lt50'
    %,..
    %'th2.6th2',...
    %{
    'th3.9th3',...
    'us3us2'...
    ,...
    'us3us2 new',...
    'th3th2 abs',...
    'us3us2 abs'
    %}
    };
%manual command line legend: legend('us3us2_new','th3th2_oldimp_abs','us3us2_oldimp_abs','Spatial LT Sigma5.2, Temporal LT Sigma4, Height0.1','Spatial LT Sigma5.2, Temporal LT Sigma50, Height0.1','Spatial TH Sigma2.6, Temporal TH Sigma2','Spatial TH Sigma3.9, Temporal TH Sigma3','Spatial US Sigma3 Alpha1.5, Temporal US Sigma2 Alpha1, Height0.5')
%'Spatial LT Sigma5.2, Temporal LT Sigma4, Height0.1',...
%'Spatial LT Sigma5.2, Temporal LT Sigma50, Height0.1',...
%'Spatial TH Sigma2.6, Temporal TH Sigma2',...
%'Spatial TH Sigma3.9, Temporal TH Sigma3',...
%'Spatial US Sigma3 Alpha1.5, Temporal US Sigma2 Alpha1, Height0.5',...
chartTitle = 'Velocity Correlations for 5 Different Paddle Correlation Kernels. Reynolds Number Constant (within 5%)';  
%chartTitle = 'Correlation Functions for Top Hat Long Tail, SpatialSigma=3.9, TemporalSigma=.3sec';

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
%workspace1 = load('statscorr_test_0810.mat');
%workspaceArray = [workspace1];
%workspaceNames = {'test'};
%chartTitle = 'Correlation Function';

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
    [ MASCc, sepvalc ] = cutoffcorrelation(workspaceArray(j).MASC, workspaceArray(j).sepval);
    
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
    hax = gca;
    ylabel('correlation   ');
    
    %loglog(sepvalc/L,MASCc,'*');
    %semilogy(sepvalc/L,MASCc);
    %semilogy(sepvalc/L,MASCc,'LineWidth',2);
    %hax = gca; 
    %ylabel('correlation   ');
    %xlabel('distance (m/L)  ');
    %xlim([0 4]);
    %ylim([0 1]);
    h = semilogy(sepvalc/workspaceArray(j).oneOverEScale,MASCc,'LineWidth',2);
    xlabel('distance (m/oneOverEScale)');
    xlim([0 4]);
    ylim([.05 1]);
    set(h, 'DisplayName', workspaceNames{j}); 
    %plot(sepvalc/L,MASCc);
    %xlabel('distance (m/L)');
    %xlim([0 4]);
    %ylim([0 1]);
    %semilogx(sepvalc/eta,MASCc);
    
    calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(j).MASvss,false);
    
    %{
    %title('Correlation Function loglog');
    title(chartTitle);

    %ncorf = fullfile(pathname, 'ncorel_trd1.5.fig');
    logcorf = fullfile(path, figurename);

    %saveas(H1, ncorf);
    saveas(H2, logcorf);
    %}
    
    % print results
    fprintf('One over e scale = %f\n',workspaceArray(j).oneOverEScale);
    fprintf('Energy dissipation rate (W) = %f \n', epsilon);
    fprintf('Kolmogorov length scale (m) = %f \n', eta);
    fprintf('Maximum fluctuation frequency (Hz) = %f \n', freq);
    hold on;
end
    
legend(workspaceNames);

rmpath(path);
rmpath(path1);
rmpath(path2);
rmpath(path3);
rmpath(path4);
rmpath(path5);
rmpath(path6);
rmpath(path7);

% play sound to alert sleeping user to end of data processing
%{
t = 0:(1/8000):0.25;
y1 = sin(2*pi*440*t);
y2 = sin(2*pi*554.37*t);
y3 = sin(2*pi*659.25*t);
y4 = sin(2*pi*880*t);
y = [y1 y2 y3 y4 y1 y1];
sound(y, 8000);
%}
