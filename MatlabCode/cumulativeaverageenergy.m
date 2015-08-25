%{
---------------------------------------------------------------------------
This program reads velocities from a .mat file and plots a cumulative
average over time. This plot represents the fluctuations of the average
over the time of data collection, and gives us an idea of the convergence
of the average velocity readings (i.e. when the curve settles down and
becomes steady, the velocity readings have converged on a stable value).

Written by Kevin Griffin, August 2015
---------------------------------------------------------------------------
%}

% load the workspace 
%workspace = load('/Users/kevin/Documents/Data/data08_17_15/statscorr_lt5.2lt50_0817_08.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_00.mat','u','MASvsm');
%workspace = load('/Users/kevin/Documents/Data/data08_14_15/statscorr_lt5.2lt50_0814_09.mat','u','MASvsm');
workspace1 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_00.mat','u','MASvsm');
workspace2 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_01.mat','u','MASvsm');
workspace3 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_02.mat','u','MASvsm');
workspace4 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_03.mat','u','MASvsm');
workspace5 = load('/Users/kevin/Documents/Data/data08_19_15/statscorr_lt5.2lt50_0819_04.mat','u','MASvsm');

workspaceArray = [workspace1, workspace2,workspace3,workspace4,workspace5];

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