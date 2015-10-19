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
sigmaArray = [2,3,4,5];

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
xlabel('Sigma');
xlim([min(sigmaArray)-1 max(sigmaArray)+1]);
ylim([0 400]);
h = gca;
set(h,'XScale','linear');
set(h,'YScale','linear');
set(h,'Fontsize', 12);
hold on;
plot(sigmaArray, turbRe);
legend('Turbulent Re');

figure();
hold on;
xlim([min(sigmaArray)-1 max(sigmaArray)+1]);
%ylim([])l
set(h,'Fontsize', 12);
title('Effect of Sigma on Length Scale and RMS Velocity at Fixed Re (within 21%)');
ylabel('Length Scales and RMS Velocity');
xlabel('Sigma');
plot(sigmaArray, taylorL);
plot(sigmaArray, LArray);
plot(sigmaArray, oneOverEScaleArray)
plot(sigmaArray, MASvssArray);
legend('Taylor Length Scale','Integral Length Scale','OneOverE Length Scale','RMS Velocity of Flow','Location','best');

rmpath(path1);
rmpath(path2);
rmpath(path3);