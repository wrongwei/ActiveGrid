pathArray = {...
    '/n/homeserver2/user3a/kevinpg/Data/data08_14_15/'... %path for 10 min lt5.2lt50
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_19_15/'... %path for 1st 40 min lt5.2lt50
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_20_15/'... %path for 2nd 40min lt5.2lt50
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_21_15/'... %path for 2nd 40min lt5.2lt50
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_25_15/'... %path for 40min lt5.2lt4 and th3.9th3
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_26_15/'... %path for 40min lt5.2lt4 and th3.9th3 and 70 min lt5.2lt50
    ,...
    '/n/homeserver2/user3a/kevinpg/Data/data08_27_15/'... %path for 2nd 40min lt5.2lt4 and th3.9th3
    };
%array of the filebases for each set of mat files
fileBaseArray = {...
    
    'energy_lt5.2lt50_0814',... %energy file for 10 min lt5.2lt50
    'energy_lt5.2lt50_0819',... %energy file for 1st 40 min lt5.2lt50
    'energy_lt5.2lt50_0820a',...%energy file for 2nd 40 min lt5.2lt50
    'energy_lt5.2lt50_0826'...  %energy file for 70 min lt5.2lt50
    
    %'energy_th3.9th3_0821a'...  %energy file for 10 min th3.9th3
    %'energy_th3.9th3_0825',...  %energy file for 40 min th3.9th3
    %'energy_th3.9th3_0827'...   %energy file for 2nd 40 min th3.9th3
    
    %'energy_lt5.2lt4_0826',...  %energy file for 40 min lt5.2lt4
    %'energy_lt5.2lt4_0827'...  %energy file for 2nd 40 min lt5.2lt4
    };
errorbars = true; 
numberOfDistances = 10; % starting with test 0
% initial parameters for curve fit (suggested: 1 3 -1)
%b0 = [1, -1];
%b0 = [0.0817, -0.9534];
% initial parameters for 3 parameter fit 
b0 = [1, 3, -1];
%b0=[1,2,-2];
samplingFrequency = 20000; % hz
paddled = 0.115;
%-------------------------------------------------------------------------
numberOfMatFiles = length(fileBaseArray);
numberOfPaths = length(pathArray);
for pathCounter = 1 : numberOfPaths
    addpath(pathArray{pathCounter});
end

%for each normalized distance, make a point on the normalized energy graph
totalLength = 0;
weightedAverageEnergyArray = zeros(numberOfDistances,1);
weightedAverageEnergyVarianceArray = zeros(numberOfDistances,1);
weightedAverageMASvsmArray = zeros(numberOfDistances,1);
weightedAverageOneOverEScaleArray = zeros(numberOfDistances,1);
%for all the mat files, take a weighted average of normalized energy
for matFileCounter = 1 : numberOfMatFiles
    filestr = strcat(fileBaseArray{matFileCounter},'.mat'); % find file name
    disp(filestr);
    vars = load(filestr,'energyArray', 'dist','rCmax','energyVarianceArray','oneOverEScaleArray','MASvsmArray');
    totalLength = totalLength + vars.rCmax;
    %weightedAverageEnergyArray;
    %vars.normEnergyArray;
    weightedAverageEnergyArray = weightedAverageEnergyArray + vars.rCmax*vars.energyArray;
    weightedAverageMASvsmArray = weightedAverageMASvsmArray + vars.rCmax*vars.MASvsmArray;
    weightedAverageOneOverEScaleArray = weightedAverageOneOverEScaleArray + vars.rCmax*vars.oneOverEScaleArray;
    if(errorbars)
        weightedAverageEnergyVarianceArray = weightedAverageEnergyVarianceArray + vars.rCmax*vars.energyVarianceArray;
    end
end
weightedAverageEnergyArray = weightedAverageEnergyArray./totalLength;
weightedAverageMASvsmArray = weightedAverageMASvsmArray./totalLength;
weightedAverageOneOverEScaleArray = weightedAverageOneOverEScaleArray./totalLength;

%plot the lengthscale versus distance
figure(2);
hold on;
plot(vars.dist-1,weightedAverageOneOverEScaleArray);
xlabel('Distance from the grid (m)');
ylabel('1/e Length Scale (m)');
%legend('lt5.2lt50','th3.9th3','lt5.2lt4');
title('Length scale vs distance from grid for the three decay kernels');

if(errorbars)
    weightedAverageEnergyVarianceArray = weightedAverageEnergyVarianceArray./totalLength;
end

if(errorbars)
    totalNumberOfCorrLengths = totalLength*weightedAverageMASvsmArray...
        ./weightedAverageOneOverEScaleArray/samplingFrequency;
    stderr = sqrt(weightedAverageEnergyVarianceArray./totalNumberOfCorrLengths);
    result = Edecfit(vars.dist-1, weightedAverageMASvsmArray, weightedAverageEnergyArray, b0, stderr);
else
    result = Edecfit(vars.dist-1, weightedAverageMASvsmArray, weightedAverageEnergyArray, b0);
end

%{
figure(2);
scatter((vars.dist-1),weightedAverageOneOverEScaleArray);
hold on;
func = @(b,x)b(1).*(x-b(2)).^b(3); %the power law y = a(x-b)^c
%initb = [1, 3, -1];
initb = b0;
[b1, R, J, covb, mse] = nlinfit((vars.dist-1)/paddled,weightedAverageOneOverEScaleArray, func, initb);
disp(b1);
real(b1);
disp(b1);
%CI = nlparci(b1, R, 'covar', covb);
%disp(CI);
x1 = linspace(0,10/paddled,10000000);
plot(x1, b1(1).*(x1-b1(2)).^b1(3));
plot(x1,paddled*ones(length(x1),1));
xlabel('Distance from Grid (m)');
ylabel('one over e scale');
set(gca,'FontSize',12);
legend('scatter of one over e scale along length of tunnel','Power Law fit a=0.4 b=1.6 c=0.1','Width of 1 paddle','location','best');
ylim([0 0.55]);
%}