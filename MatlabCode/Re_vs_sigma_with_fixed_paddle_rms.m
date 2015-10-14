%Greg's Graph #5: Effect of varying sigma on Re (with fixed Amplitude)
close all;
clear all;

% load all the workspaces you want to graph
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_12_15/');
addpath(path1);

path2 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_28_15/');
addpath(path2);

workspace1 = load('statscorr_th2.6th2_0812.mat','MASC','sepval','MASvss','oneOverEScale'); %'MASC','sepval',
workspace2 = load('statscorr_th3.9th3_rms20_0828.mat','MASC','sepval','MASvss','oneOverEScale');
workspace3 = load('statscorr_th5.2th4_rms20_0828.mat','MASC','sepval','MASvss','oneOverEScale');
workspace4 = load('statscorr_th6.5th5_rms20_0828.mat','MASC','sepval','MASvss','oneOverEScale');

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
    [ turbRe(i), taylorL(i) ] =  calculatereynoldsnumber(MASCc,sepvalc,workspaceArray(i).MASvss, true)
end

figure();
title('Effect of varying sigma on Reynolds number with constant paddle Amplitude (RMS = 20)');
ylabel('Reyonlds Number');
xlabel('Sigma');
xlim([min(sigmaArray)-1 max(sigmaArray)+1]);
ylim([0 400]);
h = gca;
set(h,'XScale','linear');
set(h,'YScale','linear');
set(h,'Fontsize', 12);
hold on;

%p = polyfit(sigmaArray,tubRe,1); %fit a second order polynomial to these 26 points. Note polyfit wants it's inputs to be row vectors, wo we use the transpose of corrVals
%x1 = linspace(0,60);
%y1 = polyval(p,x1);
%plot(x1,y1,'r')
scatter(sigmaArray,turbRe, 1000, '.');
%legend(strcat('Sigma =',sigmaArray(1)),'RMS=20','RMS=30','RMS=40','RMS=50','Location','best');

rmpath(path1);
rmpath(path2);