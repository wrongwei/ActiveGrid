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

workspaceArray = [...
    %load('statscorr_lt5.2lt50_0826_00.mat','MASC','sepval','MASvss','oneOverEScale'); % for 70min lt5.2lt50
    load('statscorr_th3.9th3_0827_08.mat','MASC','sepval','MASvss','oneOverEScale');  % for 2nd 40 min th3.9th3
    %load('statscorr_lt5.2lt4_0827_08.mat','MASC','sepval','MASvss','oneOverEScale'); % for 2nd 40min lt5.2lt4...
    ];

%manual command line legend: legend('x = 8.64 m','x = 6.23 m','x = 4.42 m','x = 3.06 m','x = 2.04 m')

% -----------------------------------------------------------------------
fprintf('Workspaces loaded in %.1f seconds.\n', round(10*toc)/10);
for j = 1 : length(workspaceArray)
    fprintf('\n-----------------------------------------------------------------------');
    fprintf('------\n Processing test \n');
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

    figure(1);
    hold on;
    set(gca, 'fontsize', 12);
    hax = gca;
    ylabel('Correlation');
    title('Correlation along the length of tunnel for th3.9th3')
    
    [turbRe, taylorL] = calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(j).MASvss,false);
    
    h = semilogy(sepvalc/workspaceArray(j).oneOverEScale,MASCc,'LineWidth',2);
    %h = semilogy(sepvalc,MASCc,'LineWidth',2);
    %h = semilogy(sepvalc*taylorL/L^3,MASCc,'LineWidth',2);
    xlabel('Distance (m/oneOverEScale)');
    %xlabel('distance (m)');
    %xlabel('distance (m*taylorL/L^3)');
    
    xlim([0 4]);
    %ylim([.05 1]);
end

rmpath(path);
rmpath(path1);
rmpath(path2);
rmpath(path3);
rmpath(path4);
rmpath(path5);
rmpath(path6);