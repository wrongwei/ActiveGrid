% Quick script to get basic statistics from energy decay tests,
% rather than running CorrCheck and loglogcorf a bajillion times
% Requirements: angle files with indices 0-9 (missing files are ok), and
% stats files generated by makeallstats with indices 0-9
% Dependencies: hwils
% Nathan Wei, August 2015

clear all;
trials = 10.0; % number of tests carried out
folder = fileparts('/Users/nathan/Documents/Data/data07_31_15/');
addpath(folder);
mean_rmsd = zeros(trials,1); % mean RMSD for a given trial
L = zeros(trials,1); % integral length scale for a given trial
Re = zeros(trials,1); % turbulence Reynolds number for a given trial
epsilon = zeros(trials,1);
nu = 1.46e-5; % kinematic viscosity of air
sum_rmsd = 0;
total_angle_files = trials;
%{
% Compute average RMSD of angles from angle files
for f = 0 : (trials - 1)
    try
        % CHANGE THIS: load angle file for test #f
        filestring1 = strcat('angles_g3g3_0731_0', num2str(f), '.txt');
        A = load(filestring1);
    catch err
        % file does not exist, so skip it and move on to the next one
        total_angle_files = total_angle_files - 1; % decrement counter
        continue;
    end
    % set up parameters
    size = size(A);
    time = size(1);
    avgangle = zeros(time,1); % stores average angle of each timeslice
    rmsd = zeros(time,1); % stores rmsd of each timeslice
    % Get mean angle of each timeslice
    for t = 1 : time
        timeslice = A(t,:);
        avgangle(t) = mean(timeslice);
    end
    % Get RMSD of each timeslice
    for t = 1 : time
        agrid = transpose(reshape(A(t,:),13,11));
        for i = 1 : 11
            for j = 1 : 13
                rmsd(t) = rmsd(t) + (agrid(i,j)-avgangle(t))^2 / 143;
            end
        end
        rmsd(t) = sqrt(rmsd(t));
    end
    mean_rmsd(f+1) = mean(rmsd);
    sum_rmsd = sum_rmsd + mean_rmsd(f+1); % get average by adding all rmsds
    fprintf('Mean RMSD of file %d = %.4f \n', f, mean_rmsd(f+1));
    clear size;
end
fprintf('Average RMSD for this run = %.4f \n', sum_rmsd/total_angle_files);
%}
% Calculating integral length scale and Reynolds number approximation
% If stats files do not exist, this code block will make them
for f = 0 : (trials - 1)
    % CHANGE THIS (should correspond to same test as angle files)
    filestring2 = strcat('statscorr_g3g3_0731_0', num2str(f), '.mat');
    load(filestring2, 'MASvss', 'MASC', 'sepval'); % std/rms, corr, ?
    fprintf('RMS velocity for file %d (m/s) = %.8f \n', f, MASvss);
    % calculate integral length scale
    L(f+1) = hwils(MASC, sepval, 2); %this is the integral length scale
    fprintf('Integral length scale for file %d (m) = %.8f \n', f, L(f+1));
    % calculate energy dissipation
    epsilon(f+1) = 0.5 * (MASvss^3) / L(f+1); % 0.5 is constant prefactor
    % calculate Kolmogorov length scale
    eta(f+1) = (nu^0.75) * (epsilon(f+1)^(-0.25));
    % calculate maximum frequency
    freq(f+1) = MASvss/eta(f+1);
    fprintf('Energy dissipation rate for file %d (W) = %.8f \n', f, epsilon(f+1));
    fprintf('Kolmogorov length scale for file %d (m) = %.8f \n', f, eta(f+1));
    fprintf('Maximum fluctuation frequency for file %d (Hz) = %.8f \n', f, freq(f+1));
    % calculate Reynolds number approximation, based on RMSD fluctuation
    % velocity (i.e. standard deviation) and integral length scale
    Re(f+1) = L(f+1)*MASvss/nu; % NEEDS TAYLOR MICROSCALE
    fprintf('Reynolds number for file %d is approximately %.4f \n\n', f, Re(f+1));
end
fprintf('Integral Length Scale range for this run: %.4f to %.4f \n', min(L), max(L));
fprintf('Reynolds Number range for this run: %.4f to %.4f\n', min(Re), max(Re));

rmpath(folder);