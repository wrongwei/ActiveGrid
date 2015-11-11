% Plot length scales vs. effective sigma for as much data as possible
% Length scales used: integral length scale, 1/e length scale
% Written by Nathan Wei and Kevin Griffin, 6 November 2015
% Dependencies: cutoffcorrelation.m, calculatereynoldsnumber.m, hwils.m,
% calculateeffectivesigma.m
close all;
paddled = 0.115; % meters
meanU = 1.45; % m/s
timeStep = 0.1; % seconds

% Easier switching between user paths
isNathan = true; % CHANGE THIS IF YOUR NAME IS KEVIN
if (isNathan)
    pathstring = '/n/homeserver2/user3a/nwei/Documents/Turbulence2015/';
else
    pathstring = '/n/homeserver2/user3a/kevinpg/Data/';
end

path = fileparts(pathstring);

% load all the workspaces you want to graph
path1 = fileparts(strcat(pathstring, 'data08_06_15/'));
addpath(path1);

path2 = fileparts(strcat(pathstring, 'data08_11_15/'));
addpath(path2);

path3 = fileparts(strcat(pathstring, 'data08_12_15/'));
addpath(path3);

path4 = fileparts(strcat(pathstring, 'data08_17_15/'));
addpath(path4);

path5 = fileparts(strcat(pathstring, 'data08_24_15/'));
addpath(path5);

path6 = fileparts(strcat(pathstring, 'data08_28_15/'));
addpath(path6);

% Start with all 29 tests in condensed summary table (add more if needed)
filenames = {...
    'lt3.9lt3_h0_0811.mat',...     1
    'lt3.9lt3_h0.05_0811.mat',...  2
    'lt3.9lt3_h0.1_0811.mat',...   3
    'lt3.9lt3_h0.2_0811.mat',...   4
    'lt3.9lt3_h0.4_0811.mat',...   5
    'lt3.9lt3_h0.8_0811.mat',...   6
    'lt6.5lt5_h0_0811.mat',...     7
    'lt6.5lt5_h0.05_0811.mat',...  8
    'lt6.5lt5_h0.1_0811.mat',...   9
    'lt6.5lt5_h0.2_0811.mat',...   10
    'lt6.5lt5_h0.4_0811.mat',...   11
    'lt6.5lt5_h0.8_0811.mat',...   12
    'th1.3th1_rms10_0806.mat',...  13
    'th1.3th1_rms20_0806.mat',...  14
    'th1.3th1_rms30_0806.mat',...  15
    'th1.3th1_rms40_0806.mat',...  16
    'th1.3th1_rms50_0806.mat',...  17
    'th2.6th2_0812.mat',...        18 
    'lt5.2lt4_0812.mat',...        19
    'lt5.2lt50_0812.mat',...       20
    'th3.9th3_0817.mat',...        21
    'th5.2th50_0824.mat',...       22
    'lt5.2lt25_0824.mat',...       23
    'lt2lt50_0824.mat',...         24
    'th5.2th4_0828.mat',...        25
    'th6.5th5_0828.mat',...        26
    'th3.9th3_rms20_0828.mat',...  27
    'th5.2th4_rms20_0828.mat',...  28
    'th6.5th5_rms20_0828.mat',...  29
    };
if (~isNathan)
    strcat('statscorr_', filenames); % should take care of statscorr vs. no statscorr issue
end

% don't load workspaces every time if they already exist in memory
loadWorkspaces = true;
if (exist('workspaceArray','var'))
    if (length(workspaceArray) == length(filenames))
        loadWorkspaces = false;
    end
end

if (loadWorkspaces)
    workspace1 = load(filenames{1},'MASC','sepval','MASvss','oneOverEScale');
    workspace2 = load(filenames{2},'MASC','sepval','MASvss','oneOverEScale');
    workspace3 = load(filenames{3},'MASC','sepval','MASvss','oneOverEScale');
    workspace4 = load(filenames{4},'MASC','sepval','MASvss','oneOverEScale');
    workspace5 = load(filenames{5},'MASC','sepval','MASvss','oneOverEScale');
    workspace6 = load(filenames{6},'MASC','sepval','MASvss','oneOverEScale');
    workspace7 = load(filenames{7},'MASC','sepval','MASvss','oneOverEScale');
    workspace8 = load(filenames{8},'MASC','sepval','MASvss','oneOverEScale');
    workspace9 = load(filenames{9},'MASC','sepval','MASvss','oneOverEScale');
    workspace10 = load(filenames{10},'MASC','sepval','MASvss','oneOverEScale');
    workspace11 = load(filenames{11},'MASC','sepval','MASvss','oneOverEScale');
    workspace12 = load(filenames{12},'MASC','sepval','MASvss','oneOverEScale');
    workspace13 = load(filenames{13},'MASC','sepval','MASvss','oneOverEScale');
    workspace14 = load(filenames{14},'MASC','sepval','MASvss','oneOverEScale');
    workspace15 = load(filenames{15},'MASC','sepval','MASvss','oneOverEScale');
    workspace16 = load(filenames{16},'MASC','sepval','MASvss','oneOverEScale');
    workspace17 = load(filenames{17},'MASC','sepval','MASvss','oneOverEScale');
    workspace18 = load(filenames{18},'MASC','sepval','MASvss','oneOverEScale');
    workspace19 = load(filenames{19},'MASC','sepval','MASvss','oneOverEScale');
    workspace20 = load(filenames{20},'MASC','sepval','MASvss','oneOverEScale');
    workspace21 = load(filenames{21},'MASC','sepval','MASvss','oneOverEScale');
    workspace22 = load(filenames{22},'MASC','sepval','MASvss','oneOverEScale');
    workspace23 = load(filenames{23},'MASC','sepval','MASvss','oneOverEScale');
    workspace24 = load(filenames{24},'MASC','sepval','MASvss','oneOverEScale');
    workspace25 = load(filenames{25},'MASC','sepval','MASvss','oneOverEScale');
    workspace26 = load(filenames{26},'MASC','sepval','MASvss','oneOverEScale');
    workspace27 = load(filenames{27},'MASC','sepval','MASvss','oneOverEScale');
    workspace28 = load(filenames{28},'MASC','sepval','MASvss','oneOverEScale');
    workspace29 = load(filenames{29},'MASC','sepval','MASvss','oneOverEScale');

    workspaceArray = [workspace1,workspace2,workspace3,workspace4,workspace5,...
        workspace6,workspace7,workspace8,workspace9,workspace10,workspace11...
        ,workspace12,workspace13,workspace14,workspace15,workspace16,workspace17...
        ,workspace18,workspace19,workspace20,workspace21,workspace22,workspace23...
        ,workspace24,workspace25,workspace26,workspace27,workspace28,workspace29]; 

end
disp('Workspaces loaded!');

% String processing hacks
filenames{24} = strrep(filenames{24},'lt2','lt2.0');
filenames = strrep(filenames, 'statscorr_', ''); % take off statscorrs if necessary
filenames{1} = strrep(filenames{1},'h0','h0.0'); % height finding code needs #.# format
filenames{7} = strrep(filenames{7},'h0','h0.0');

% set up necessary arrays
turbRe = zeros(length(workspaceArray), 1);
taylorL = zeros(length(workspaceArray), 1);
LArray = zeros(length(workspaceArray), 1);
MASvssArray = zeros(length(workspaceArray), 1);
oneOverEScaleArray = zeros(length(workspaceArray), 1);
spatialSigmas = zeros(length(workspaceArray), 1);
temporalSigmas = zeros(length(workspaceArray), 1);
effectiveSigmas = zeros(length(workspaceArray), 1);

for i = 1 : length(workspaceArray)
    % use helper functions to calculate turbulent Reynolds number and
    % length scales
    disp(strcat('Processing test ', num2str(i)));
    [ MASCc, sepvalc ] = cutoffcorrelation(workspaceArray(i).MASC, workspaceArray(i).sepval);
    [ turbRe(i), taylorL(i) ] =  calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(i).MASvss, false);
    LArray(i) = hwils(MASCc,sepvalc,2);
    MASvssArray(i) = workspaceArray(i).MASvss;
    oneOverEScaleArray(i) = workspaceArray(i).oneOverEScale;
    [spatialSigmas(i), temporalSigmas(i), effectiveSigmas(i)] =...
        calculateeffectivesigma(filenames{i}, paddled, meanU, timeStep);
end

%disp(oneOverEScaleArray);
%disp(LArray);

% plot length scales with effective sigma
% Plot length scale vs. effective sigma
H1 = figure(1);
set(gca, 'fontsize', 12);
hax = gca;
scatter(effectiveSigmas, LArray, 'o');
hold on;
scatter(effectiveSigmas, oneOverEScaleArray, '*');
scatter(effectiveSigmas, taylorL, '.');
title('Length Scales vs. Effective Sigma');
xlabel('Effective Sigma (m)');
ylabel('Length Scales (m)');
xlim([0 30]);
ylim([0 4]);
set(hax,'XScale','log');
legend('Integral Length Scale','1/e Length Scale','Taylor Length Scale',...
    'location','northwest');
%text(effectiveSigmas, LArray, filenames);
%text(effectiveSigmas, oneOverEScaleArray, filenames);
%text(effectiveSigmas, taylorL, filenames);
plot1 = fullfile(path, 'length_scales_vs_effective_sigma.fig');
saveas(H1, plot1);

% clean up
rmpath(path1);
rmpath(path2);
rmpath(path3);
rmpath(path4);
rmpath(path5);
rmpath(path6);