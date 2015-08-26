% Quick script to get actual paddle RMS angles for multiple files,
% rather than running CorrCheck a bajillion times.
% Requirements: angle files in a directory (specify the search string using
%   regular expressions within the dir() function argument
% Dependencies: none
% Nathan Wei, August 2015

clear all;
folder = fileparts('/Users/nathan/Documents/Data/data08_21_15/');
addpath(folder);

% find all angle files in a given folder (angle files start with "angles_")
files = dir('/Users/nathan/Documents/Data/data08_21_15/angles_*');

% Compute average RMSD of angles from angle files
for f = 1 : length(files)
    A = load(files(f).name);
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
    fprintf('Mean RMSD of file %s = %f \n', files(f).name, mean(rmsd));
    clear size;
end

rmpath(folder);