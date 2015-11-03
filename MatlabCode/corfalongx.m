% plots the normalized (by integral length scale) correlation function at
% different x postions fr the three kernels used in decay.

% Where is the path?-----------------------------------------------------
path = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_14_15/');
addpath(path);
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_19_15/');
addpath(path1);
path2 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_20_15/');
addpath(path2);
path3 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_21_15/');
addpath(path3);
path4 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_25_15/');
addpath(path4);
path5 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_26_15/');
addpath(path5);
path6 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_27_15/');
addpath(path6);
% load all the workspaces you want to graph. Put each one in a varaible,
% and then put all of those variables into the array below named
% workspaceArray
% -----------------------------------------------------------------------

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('*****************************************************************************\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Loading workspaces... ');
tic;

%workspace1 = load('statscorr_lt5.2lt50_0814.mat','MASC','sepval','MASvss','oneOverEScale'); % for 10 min lt5.2lt50
%workspace2 = load('statscorr_lt5.2lt50_0819.mat','MASC','sepval','MASvss','oneOverEScale'); % for 40min lt5.2lt50
%workspace3 = load('statscorr_lt5.2lt50_0820a.mat','MASC','sepval','MASvss','oneOverEScale'); % for 2nd 40min lt5.2lt50
workspace4 = load('statscorr_lt5.2lt50_0826_08.mat','MASC','sepval','MASvss','oneOverEScale'); % for 70min lt5.2lt50
%{
workspace1 = load('statscorr_th3.9th3_0821a.mat','MASC','sepval','MASvss','oneOverEScale'); % for 10 min th3.9th3
workspace2 = load('statscorr_th3.9th3_0825.mat','MASC','sepval','MASvss','oneOverEScale');  % for 40min th3.9th3 
workspace3 = load('statscorr_th3.9th3_0827.mat','MASC','sepval','MASvss','oneOverEScale');  % for 2nd 40 min th3.9th3
%}
%{
workspace1 = load('statscorr_lt5.2lt4_0826.mat','MASC','sepval','MASvss','oneOverEScale'); % for 40min lt5.2lt4
workspace2 = load('statscorr_lt5.2lt4_0827.mat','MASC','sepval','MASvss','oneOverEScale'); % for 2nd 40min lt5.2lt4
%}

workspaceArray = [...
    %workspace1...
    %,...
    %workspace2...
    %,...
    %workspace3...
    %,...
    workspace4...
    ];

workspaceNames = {...
    %'lt5.2lt50'...
    %,...
    'th3.9th3'...
    %,..
    %'lt5.lt4'...
    };
%manual command line legend: legend('us3us2_new','th3th2_oldimp_abs','us3us2_oldimp_abs','Spatial LT Sigma5.2, Temporal LT Sigma4, Height0.1','Spatial LT Sigma5.2, Temporal LT Sigma50, Height0.1','Spatial TH Sigma2.6, Temporal TH Sigma2','Spatial TH Sigma3.9, Temporal TH Sigma3','Spatial US Sigma3 Alpha1.5, Temporal US Sigma2 Alpha1, Height0.5')
%'Spatial LT Sigma5.2, Temporal LT Sigma4, Height0.1',...
%'Spatial LT Sigma5.2, Temporal LT Sigma50, Height0.1',...
%'Spatial TH Sigma2.6, Temporal TH Sigma2',...
%'Spatial TH Sigma3.9, Temporal TH Sigma3',...
%'Spatial US Sigma3 Alpha1.5, Temporal US Sigma2 Alpha1, Height0.5',...
chartTitle = 'Velocity Correlations for 5 Different Paddle Correlation Kernels. Reynolds Number Constant (within 5%)';  
%chartTitle = 'Correlation Functions for Top Hat Long Tail, SpatialSigma=3.9, TemporalSigma=.3sec';

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
    
    [turbRe, taylorL] = calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(j).MASvss,false);
    
    h = semilogy(sepvalc/workspaceArray(j).oneOverEScale,MASCc,'LineWidth',2);
    %h = semilogy(sepvalc,MASCc,'LineWidth',2);
    %h = semilogy(sepvalc*taylorL/L^3,MASCc,'LineWidth',2);
    xlabel('distance (m/oneOverEScale)');
    %xlabel('distance (m)');
    %xlabel('distance (m*taylorL/L^3)');
    
    xlim([0 4]);
    ylim([.05 1]);
    set(h, 'DisplayName', workspaceNames{j}); 
    
    %plot(sepvalc/L,MASCc);
    %xlabel('distance (m/L)');
    %xlim([0 4]);
    %ylim([0 1]);
    %semilogx(sepvalc/eta,MASCc);
    
   
    
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