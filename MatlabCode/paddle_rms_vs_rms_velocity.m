%Greg's Graph #4: Effect of varying paddle RMS (amplitude) on RMS Velocity
clear all;

% load all the workspaces you want to graph
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_06_15/');
path2 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_07_15/');
path3 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_11_15/old_statscorr/');
path4 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_28_15/');
addpath(path1);
addpath(path2);
addpath(path3);
addpath(path4);

% enter the arrays of rms angles into this section
angleArray1 = [9.808886;19.505719;28.51177;35.983436;41.412413];
angleArray2 = [16; 17.5; 20];
angleArray3 = [14; 17; 20];
angleArray4 = [40; 38; 35; 27; 19; 17];
angleArray5 = [30; 21; 25; 27];
angleArray6 = [10; 40];
angleArray7 = [0];

paddleAmplitude = {angleArray1,angleArray2,angleArray3,angleArray4,angleArray5,angleArray6,angleArray7};

workspace1 = load('statscorr_th1.3th1_rms10_tr0.1_0806.mat','MASvss','MASvsm');
workspace2 = load('statscorr_th1.3th1_rms20_tr0.1_0806.mat','MASvss','MASvsm');
workspace3 = load('statscorr_th1.3th1_rms30_tr0.1_0806.mat','MASvss','MASvsm');
workspace4 = load('statscorr_th1.3th1_rms40_tr0.1_0806.mat','MASvss','MASvsm');
workspace5 = load('statscorr_th1.3th1_rms50_tr0.1_0806.mat','MASvss','MASvsm');
workspaceArray1 = [workspace1, workspace2,  workspace3, workspace4, workspace5]; 

workspace1 = load('statscorr_th5.2th4_a_0828.mat','MASvss','MASvsm');
workspace2 = load('statscorr_th5.2th4_b_0828.mat','MASvss','MASvsm');
workspace3 = load('statscorr_th5.2th4_rms20_0828.mat','MASvss','MASvsm');
workspaceArray2 = [workspace1, workspace2, workspace3];

workspace1 = load('statscorr_th6.5th5_a_0828.mat','MASvss','MASvsm');
workspace2 = load('statscorr_th6.5th5_b_0828.mat','MASvss','MASvsm');
workspace3 = load('statscorr_th6.5th5_rms20_0828.mat','MASvss','MASvsm');
workspaceArray3 = [workspace1, workspace2, workspace3];

workspace1 = load('statscorr_lt1.3lt1_h0_0807.mat','MASvss','MASvsm');
workspace2 = load('statscorr_lt1.3lt1_h0.05_0807.mat','MASvss','MASvsm');
workspace3 = load('statscorr_lt1.3lt1_h0.1_0807.mat','MASvss','MASvsm');
workspace4 = load('statscorr_lt1.3lt1_h0.2_0807.mat','MASvss','MASvsm');
workspace5 = load('statscorr_lt1.3lt1_h0.4_0807.mat','MASvss','MASvsm');
workspace6 = load('statscorr_lt1.3lt1_h0.8_0807.mat','MASvss','MASvsm');
workspaceArray4 = [workspace1, workspace2,  workspace3, workspace4, workspace5,workspace6];

workspace1 = load('statscorr_lt5.2lt30_1_0811.mat','MASvss','MASvsm');
workspace2 = load('statscorr_lt5.2lt30_2_0811.mat','MASvss','MASvsm');
workspace3 = load('statscorr_lt5.2lt30_3_0811.mat','MASvss','MASvsm');
workspace4 = load('statscorr_lt5.2lt30_4_0811.mat','MASvss','MASvsm');
workspaceArray5 = [workspace1, workspace2,  workspace3, workspace4];

workspace1 = load('statscorr_tr5.2tr50_rms10_0828.mat','MASvss','MASvsm');
workspace2 = load('statscorr_tr5.2tr50_rms40_0828.mat','MASvss','MASvsm');
workspaceArray6 = [workspace1, workspace2];

workspace1 = load('statscorr_opengrid_0828.mat','MASvss','MASvsm');
workspaceArray7 = [workspace1];

datasets = {workspaceArray1,workspaceArray2,workspaceArray3,workspaceArray4,workspaceArray5,workspaceArray6,workspaceArray7};
%Call the taylor Re function (returns taylor length scale and Re_turb)

figure(1);
title('RMS Velocity Fluctuations vs Paddle Amplitude (RMS)');
ylabel('RMS Fluctuation Velocity / Mean Velocity');
xlabel('Paddle RMS Tip Speed / Mean Velocity');
xlim([0 0.4]);
ylim([0 0.065]);
h = gca;
set(h,'XScale','linear');
set(h,'YScale','linear');
set(h,'Fontsize', 12);
hold on;
for j = 1:length(datasets)
    colorChoice = rand(1,3); %choose plotting color
    clear meanVelocityArray rmsVelocityArray;
    rmsVelocityArray = zeros(length(paddleAmplitude{j}),1);
    meanVelocityArray = zeros(length(paddleAmplitude{j}),1);
    for i = 1:length(datasets{j})
        rmsVelocityArray(i) = datasets{j}(i).MASvss;
        meanVelocityArray(i,1) = datasets{j}(i).MASvsm;
    end

    %scatter(paddleAmplitude,rmsVelocityArray, 1000, '.r');
    meanVelocityArray
    padTipLength = 0.115*sind(45); %length to tip of paddle
    updateFreq = 0.1; % number of grids per second (1/s)
    plot(paddleAmplitude{j}*padTipLength*updateFreq./meanVelocityArray,rmsVelocityArray./meanVelocityArray,'.','MarkerSize', 25,'color',colorChoice);
    %p = polyfit(paddleAmplitude,rmsVelocityArray,1);
    p = polyfit(paddleAmplitude{j}*padTipLength*updateFreq./meanVelocityArray,rmsVelocityArray./meanVelocityArray,1);
    x1 = linspace(0,60);
    y1 = polyval(p,x1);
    plot(x1,y1,'color',colorChoice)
end
legend('','th1.3th1','','th5.2th4','','th6.5th5','lt1.3lt1 with','h not constant','','lt5.2lt30','','tr5.2tr50','','open grid','Location','best');
rmpath(path1);
rmpath(path2);
rmpath(path3);
rmpath(path4);