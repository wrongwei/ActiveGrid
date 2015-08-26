%------------------------------------------------------------------------------
% MODIFY THIS SECTION
%------------------------------------------------------------------------------
%load the workspaces
workspace = load('/Users/kevin/Documents/Data/data08_24_15/statscorr_th5.2th50_0824.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_00.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_09.mat','u','MASvsm');
%workspace1 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_00.mat','u','MASvsm');
%workspace2 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_01.mat','u','MASvsm');
%workspace3 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_02.mat','u','MASvsm');
%workspace4 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_03.mat','u','MASvsm');
%workspace5 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_04.mat','u','MASvsm');

workspaceArray = [workspace];%, workspace2,workspace3,workspace4,workspace5];

legendNames = {'lt5.2lt50 00'};%;'lt5.2lt50 01';'lt5.2lt50 02';'lt5.2lt50 03';'lt5.2lt50 04'};
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
for i = 1 : length(workspaceArray)
    samplesPerGroup = 200000;
    samplesPerMinute = 20000*60;
    cumulativeAverage = zeros(length(workspaceArray(i).u)/samplesPerGroup,1);
    %movingAverageArray = zeros(length(workspaceArray(i).u)/samplesPerGroup,1);
    numberOfMinutes = (1:length(workspaceArray(i).u)/samplesPerGroup);
    numberOfMinutes = numberOfMinutes*samplesPerGroup/samplesPerMinute;

    %{
    for j = 1 : length(workspaceArray(i).u)/samplesPerGroup
        cumulativeAverage(j) =  mean((workspaceArray(i).u(1:j*samplesPerGroup)-workspaceArray(i).MASvsm).^2); 
    end
    %}


    %special case for j=1
    cumulativeAverage(1) = mean((workspaceArray(i).u(1:samplesPerGroup)-workspaceArray(i).MASvsm).^2);
    fprintf('cumulativeAverage(1) = %f\n',cumulativeAverage(1));
    movingAverage = 0;
    %movingAverageArray(1) = cumulativeAverage(1);
    %calculate the new average energy
    for j = 2 : length(workspaceArray(i).u)/samplesPerGroup
        movingAverage = mean((workspaceArray(i).u((j-1)*samplesPerGroup:j*samplesPerGroup)-workspaceArray(i).MASvsm).^2);
        %movingAverageArray(j) = movingAverage;
        cumulativeAverage(j) = (cumulativeAverage(j-1)*(j-1) + movingAverage) / (j);
    end


    %plot the cumulative average
    figure(6);
    plot(numberOfMinutes, cumulativeAverage/workspaceArray(i).MASvsm^2);
    xlabel('Minutes of Data');
    ylabel('Normalized Energy');
    %ylim([10^-3 3.5*10^-3]);
    title('Convergence of cumulative average energy over 10 minutes');
    hold on;
    %{
    %plot the moving average
    figure(5);
    plot(numberOfMinutes, movingAverageArray/workspaceArray(i).MASvsm^2);
    xlabel('Minutes of Data');
    ylabel('Normalized Energy');
    title('Moving average energy over 10 minutes');
    %}
end
legend(legendNames);