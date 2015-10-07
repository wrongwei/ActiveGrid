%Greg's Graph #4: Effect of varying paddle RMS (amplitude) on RMS Velocity

paddleAmplitude = [9.808886;19.505719;28.51177;35.983436;41.412413];

% load all the workspaces you want to graph
path1 = fileparts('/n/homeserver2/user3a/kevinpg/Data/data08_06_15/');
addpath(path1);

workspace1 = load('statscorr_th1.3th1_rms10_tr0.1_0806.mat','MASvss','oneOverEScale'); %'MASC','sepval',
workspace2 = load('statscorr_th1.3th1_rms20_tr0.1_0806.mat','MASvss','oneOverEScale');
workspace3 = load('statscorr_th1.3th1_rms30_tr0.1_0806.mat','MASvss','oneOverEScale');
workspace4 = load('statscorr_th1.3th1_rms40_tr0.1_0806.mat','MASvss','oneOverEScale');
workspace5 = load('statscorr_th1.3th1_rms50_tr0.1_0806.mat','MASvss','oneOverEScale');

workspaceArray = [...
    workspace1...
    ,...
    workspace2...
    ,...
    workspace3...
    ,...
    workspace4...
    ,...
    workspace5...
    ];   

%Call the taylor Re function (returns taylor length scale and Re_turb)

figure(1);
title('Effect of varying paddle RMS (amplitude) on RMS Velocity');
ylabel('RMS Velocity');
xlabel('Paddle Amplitude (RMS)');
xlim([0 60]);
%ylim([]);
h = gca;
set(h,'XScale','linear');
set(h,'YScale','linear');
set(h,'Fontsize', 14);
hold on;
rmsVelocityArray = zeros(length(paddleAmplitude),1);
for i = 1:length(workspaceArray)
    rmsVelocityArray(i) = workspaceArray(i).MASvss;
end

p = polyfit(paddleAmplitude,rmsVelocityArray,1); %fit a second order polynomial to these 26 points. Note polyfit wanted both vectors to be row vectors, I transpose corrVals
x1 = linspace(0,60);
y1 = polyval(p,x1);
plot(x1,y1,'r')
scatter(paddleAmplitude,rmsVelocityArray, 1000, '.');
legend('RMS=10','RMS=20','RMS=30','RMS=40','RMS=50','Location','best');

rmpath(path1);