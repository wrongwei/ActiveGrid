% load the workspace 
%workspace = load('/Users/kevin/Documents/Data/data08_17_15/statscorr_lt5.2lt50_0817_08.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_00.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_09.mat','u','MASvsm');
workspace = load('/Users/kevin/Documents/Data/data08_18_15/statscorr_lt5.2lt50_0818_05.mat','u','MASvsm');
samplesPerGroup = 200000;
samplesPerMinute = 20000*60;
cumulativeAverage = zeros(length(workspace.u)/samplesPerGroup,1);
%movingAverageArray = zeros(length(workspace.u)/samplesPerGroup,1);
numberOfMinutes = (1:length(workspace.u)/samplesPerGroup);
numberOfMinutes = numberOfMinutes*samplesPerGroup/samplesPerMinute;

%{
for j = 1 : length(workspace.u)/samplesPerGroup
    cumulativeAverage(j) =  mean((workspace.u(1:j*samplesPerGroup)-workspace.MASvsm).^2); 
end
%}

%special case for j=1
cumulativeAverage(1) = mean((workspace.u(1:samplesPerGroup)-workspace.MASvsm).^2);
fprintf('cumulativeAverage(1) = %f\n',cumulativeAverage(1));
movingAverage = 0;
%movingAverageArray(1) = cumulativeAverage(1);
%calculate the new average energy
for j = 2 : length(workspace.u)/samplesPerGroup
    movingAverage = mean((workspace.u((j-1)*samplesPerGroup:j*samplesPerGroup)-workspace.MASvsm).^2);
    %movingAverageArray(j) = movingAverage;
    cumulativeAverage(j) = (cumulativeAverage(j-1)*(j-1) + movingAverage) / (j);
end

%plot the cumulative average
figure(6);
plot(numberOfMinutes, cumulativeAverage/workspace.MASvsm^2);
xlabel('Minutes of Data');
ylabel('Normalized Energy');
title('Convergence of cumulative average energy over 10 minutes');
%{
%plot the moving average
figure(5);
plot(numberOfMinutes, movingAverageArray/workspace.MASvsm^2);
xlabel('Minutes of Data');
ylabel('Normalized Energy');
title('Moving average energy over 10 minutes');
%}