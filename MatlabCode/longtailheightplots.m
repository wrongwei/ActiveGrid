% Generates plots of long tail data with varying tail heights
% Raw data from 11 August 2015
% Dependencies: hwils.m
% Written by Kevin Griffin and Nathan Wei, October 2015

% Specify path to data folder ---------------------------------------------
path = fileparts('/n/homeserver2/user3a/nwei/Documents/Turbulence2015/data08_11_15/');
addpath(path);

% load all the workspaces you want to graph. Put each one in a varaible,
% and then put all of those variables into the array below named
% workspaceArray
% -----------------------------------------------------------------------

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('*****************************************************************************\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Loading workspaces... ');
tic;

workspace01 = load('lt3.9lt3_h0_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace02 = load('lt3.9lt3_h0.05_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace03 = load('lt3.9lt3_h0.1_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace04 = load('lt3.9lt3_h0.2_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace05 = load('lt3.9lt3_h0.4_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace06 = load('lt3.9lt3_h0.8_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace07 = load('lt6.5lt5_h0_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace08 = load('lt6.5lt5_h0.05_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace09 = load('lt6.5lt5_h0.1_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace10 = load('lt6.5lt5_h0.2_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace11 = load('lt6.5lt5_h0.4_0811.mat','MASC','sepval','MASvss','oneOverEScale');
workspace12 = load('lt6.5lt5_h0.8_0811.mat','MASC','sepval','MASvss','oneOverEScale');

workspaceArray = [workspace01, workspace02, workspace03, workspace04, ...
    workspace05, workspace06, workspace07, workspace08, workspace09, ...
    workspace10, workspace11, workspace12];

workspaceNames = {'lt3.9lt3 h0.0', 'lt3.9lt3 h0.05', 'lt3.9lt3 h0.1', 'lt3.9lt3 h0.2', ...
    'lt3.9lt3 h0.4', 'lt3.9lt3 h0.8', 'lt6.5lt5 h0.0', 'lt3.9lt3 h0.05', ...
    'lt6.5lt5 h0.1', 'lt6.5lt5 h0.2', 'lt6.5lt5 h0.4', 'lt6.5lt5 h0.8', };

% preallocating needed arrays for speed
L = zeros(length(workspaceArray), 1);
turbRe = zeros(length(workspaceArray), 1);
taylorL = zeros(length(workspaceArray), 1);

% This change in involved prefixed the loaded workspace variables with workspaceArray(j).
% These variables include MASC, MASvss, sepval, etc.
% -----------------------------------------------------------------------
if (length(workspaceNames) ~= length(workspaceArray)) % validation
    error('Array size mismatch between workspaceNames and workspaceArray');
end

fprintf('Workspaces loaded in %.1f seconds.\n', round(10*toc)/10);

H1 = figure(1);

for j = 1 : length(workspaceArray)
    fprintf('\n-----------------------------------------------------------------------');
    fprintf('------\n Processing test %s \n',workspaceNames{j});
    
    L(j) = hwils(workspaceArray(j).MASC,workspaceArray(j).sepval,2); %this is the integral length scale

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
    epsilon = 0.5 * (workspaceArray(j).MASvss^3) / L(j); % 0.5 is constant prefactor
    % calculate Kolmogorov length scale
    eta = (nu^0.75) * (epsilon^(-0.25));
    % calculate maximum frequency
    freq = workspaceArray(j).MASvss/eta;
    
    % Plot correlation function on normal/loglog axes. with x-axis normalized
    % by integral lengthscale
    % gives a vertical line on where the integral length scale is on the
    % correlation function.

    samplePoints = 5:30; % make a vector of the first 26 points excluding noise
    samplePoints2 = sepvalc(samplePoints)/L(j);
    corrVals = MASCc(samplePoints); % correlation value of the 26 sample points
    p = polyfit(samplePoints2,corrVals',2); %fit a second order polynomial to these 26 points. Note polyfit wanted both vectors to be row vectors, I transpose corrVals
    % print standard deviation = RMS velocity
    fprintf('RMS velocity (m/s) = %f \n', workspaceArray(j).MASvss);
    fprintf('Integral Length Scale = %f\n', L(j));
    % the x-intercept of polyfit p is the taylor length scale
    taylorL(j) = max(roots(p))*L(j);
    fprintf('Taylor length scale = %f\n', taylorL(j));
    %fprintf('Comment: I multipled by L because of the scaling of the graph\n');
    turbRe(j) = workspaceArray(j).MASvss*taylorL(j)/nu;
    fprintf('Turbulent Reynolds Number = %f\n', turbRe(j));
    
    % print results
    fprintf('One over e scale = %f\n',workspaceArray(j).oneOverEScale);
    fprintf('Energy dissipation rate (W) = %f \n', epsilon);
    fprintf('Kolmogorov length scale (m) = %f \n', eta);
    fprintf('Maximum fluctuation frequency (Hz) = %f \n', freq);
    hold on;
    
    % plot correlation functions
    
    hold on;

    set(gca, 'fontsize', 12);
    hax = gca;
    ylabel('Correlation');

    h = semilogy(sepvalc/workspaceArray(j).oneOverEScale,MASCc,'LineWidth',2);
    xlabel('distance (m/oneOverEScale)');
    xlim([0 4]);
    ylim([0 1]);
    set(h, 'DisplayName', workspaceNames{j}); 
end

% Finish plotting correlation functions
legend(workspaceNames);
title('Correlation Functions for Various Long Tail Heights');
plot1 = fullfile(path, 'corrfs_ltheights.fig');
saveas(H1, plot1);

% Plot integral or 1/e length scale vs. long tail height
heights = [0.0 0.05 0.1 0.2 0.4 0.8];
lt3_9lt3_lengthscales = zeros(length(workspaceArray)/2, 1);
lt6_5lt5_lengthscales = zeros(length(workspaceArray)/2, 1);
for i = 1 : length(workspaceArray)/2
    lt3_9lt3_lengthscales(i) = workspaceArray(i).oneOverEScale;
end
for i = 1 : length(workspaceArray)/2
    lt6_5lt5_lengthscales(i) = workspaceArray(i+6).oneOverEScale;
end
H2 = figure(2);
set(gca, 'fontsize', 12);
hax = gca;
plot(heights, lt3_9lt3_lengthscales);
hold on;
plot(heights, lt6_5lt5_lengthscales);
title('1/e Length Scale vs. Long Tail Height');
xlabel('Long Tail Height');
ylabel('1/e Length Scale');
xlim('auto');
ylim('auto');
legend('lt3.9lt3','lt6.5lt5');
plot2 = fullfile(path, '1overE_vs_ltheights.fig');
saveas(H2, plot2);

% Plot length scale vs. effective sigma
sigma_eff_1 = zeros(length(heights), 1);
sigma_eff_2 = zeros(length(heights), 1);
for i = 1 : length(heights)
    for x = -7 : 7
        for y = -7 : 7
            dist = sqrt(x^2 + y^2);
            if (dist <= 3.9 && dist ~= 0)
                sigma_eff_1(i) = sigma_eff_1(i) + heights(i);
            end
            if (dist == 0)
                sigma_eff_1(i) = sigma_eff_1(i) + 1;
            end
        end
    end
end
for i = 1 : length(heights)
    for x = -7 : 7
        for y = -7 : 7
            dist = sqrt(x^2 + y^2);
            if (dist <= 6.5 && dist ~= 0)
                sigma_eff_2(i) = sigma_eff_2(i) + heights(i);
            end
            if (dist == 0)
                sigma_eff_2(i) = sigma_eff_2(i) + 1;
            end
        end
    end
end
H3 = figure(3);
set(gca, 'fontsize', 12);
hax = gca;
plot(sigma_eff_1, lt3_9lt3_lengthscales);
hold on;
plot(sigma_eff_2, lt6_5lt5_lengthscales);
title('1/e Length Scale vs. Long Tail Height');
xlabel('Effective Sigma');
ylabel('1/e Length Scale');
xlim('auto');
ylim('auto');
legend('lt3.9lt3','lt6.5lt5');
plot3 = fullfile(path, '1overE_vs_sigma_eff.fig');
saveas(H3, plot3);
    
%legend(workspaceNames);
rmpath(path);
