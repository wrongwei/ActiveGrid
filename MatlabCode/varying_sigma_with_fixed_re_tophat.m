%Greg's Graph #2?: Effect of varying sigma on Lengths Scale? (with fixed Re)
close all;
clear all;

% load all the workspaces you want to graph
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_12_15/');
addpath(path1);

path2 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_17_15/');
addpath(path2);

path3 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_28_15/');
addpath(path3);

workspace1 = load('statscorr_th2.6th2_0812.mat','MASC','sepval','MASvss','oneOverEScale'); %'MASC','sepval',
workspace2 = load('statscorr_th3.9th3_0817.mat','MASC','sepval','MASvss','oneOverEScale');
workspace3 = load('statscorr_th5.2th4_b_0828.mat','MASC','sepval','MASvss','oneOverEScale');
workspace4 = load('statscorr_th6.5th5_b_0828.mat','MASC','sepval','MASvss','oneOverEScale');

workspaceArray = [...
    workspace1...
    ,...
    workspace2...
    ,...
    workspace3...
    ,...
    workspace4...
    ];   

% will need touse the relevant sigma factor
spatialSigmaArray = [2.6,3.9,5.2,6.5];
temporalSigmaArray = [2,3,4,5];
for i = 1 : length(spatialSigmaArray)
    effectiveSigmaArray(i) = (8 * spatialSigmaArray(i)^2*temporalSigmaArray(i))^(1/3);
end

%rmsVelocityArray = zeros(length(paddleAmplitude),1);
for i = 1:length(workspaceArray)
    %rmsVelocityArray(i) = workspaceArray(i).MASvss;
    [ MASCc, sepvalc ] = cutoffcorrelation(workspaceArray(i).MASC, workspaceArray(i).sepval);
    [ turbRe(i), taylorL(i) ] =  calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(i).MASvss, false);
    turbRe %print
    LArray(i) = hwils(MASCc,sepvalc,2);
    MASvssArray(i) = workspaceArray(i).MASvss;
    oneOverEScaleArray(i) = workspaceArray(i).oneOverEScale;
end

figure();
title('Attempts to Fix Re by Changing Paddle Amplitude');
ylabel('Reyonlds Number');
xlabel('Effective Sigma');
xlim([min(effectiveSigmaArray)-1 max(effectiveSigmaArray)+1]);
ylim([0 400]);
h = gca;
set(h,'XScale','linear');
set(h,'YScale','linear');
set(h,'Fontsize', 12);
hold on;
plot(effectiveSigmaArray, turbRe);
legend('Turbulent Re');

figure();
hold on;
xlim([min(effectiveSigmaArray)-1 max(effectiveSigmaArray)+1]);
%ylim([])l
set(h,'Fontsize', 12);
title('Length Scale vs effective sigma for tophats of varying widths (Re constant within 21%)');
ylabel('Length Scales');
xlabel('Effective Sigma');
%plot(effectiveSigmaArray, taylorL);
plot(effectiveSigmaArray, LArray);
plot(effectiveSigmaArray, oneOverEScaleArray)
legend('Integral Length Scale','OneOverE Length Scale','Location','best');

figure();
hold on;
xlim([min(effectiveSigmaArray)-1 max(effectiveSigmaArray)+1]);
ylim([0 0.07]);
set(h,'Fontsize', 12);
title('RMS Velocity vs effective sigma for tophats of varying widths (Re constant within 21%)');
ylabel('RMS Velocity');
xlabel('Effective Sigma');
plot(effectiveSigmaArray, MASvssArray);
legend('RMS Velocity of Flow','Location','best');

rmpath(path1);
rmpath(path2);
rmpath(path3);