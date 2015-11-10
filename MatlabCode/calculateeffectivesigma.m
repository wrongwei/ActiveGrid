% Calculate effective sigma (cube root of volume of kernel) along lines of
% normalization calculations in runcorr method (algo.cpp and algo3d.cpp)
% Takes filename, paddle length, average free stream velocity, and time
% scale of grid as argument, extracts sigmas from it, and calculates
% effective sigma (returns all 3 parameters)
% Written by Nathan Wei and Kevin Griffin, 10 November 2015
% Dependencies: none
function [spatialS, temporalS, effectiveS] = ...
    calculateeffectivesigma(filename, paddled, meanU, timeStep)
    testID = strsplit(filename, '_'); % remove test ID from filename
    testID = testID{1}; % thank goodness for matlab lazy typing!
    spatialS = str2double(testID(3:5)); % assuming format xx#.#yy#####
    temporalS = str2double(testID(8:end));
    kernelID = testID(1:2); % assuming xx = yy (same kernel type in both dimensions)
    isLT = false;
    isTH = false;
    if (strcmp(kernelID, 'lt')) % if long tail kernel
        % isolate height
        isLT = true;
        heightString = filename(strfind(filename, 'h') : strfind(filename, '.mat'));
        heightString = strsplit(heightString, '_');
        heightString = heightString{1};
        height = str2double(heightString(2:end));
        if (isnan(height))
            height = 0.1; % unspecified heights are 0.1
        end
    elseif (strcmp(kernelID, 'th')) % if top hat kernel
        isTH = true;
    else
        effectiveS = -99; % something went wrong; we should only have lt and th
        return;
    end
    effectiveS = 0;
    boundS = ceil(spatialS);
    boundT = ceil(temporalS);
    % calculate effective sigma as in algo3d.cpp
    for i = -boundS : boundS
        for j = -boundS : boundS
            for t = -boundT : boundT
                dist = sqrt(i^2 + j^2);
                seed = 1;
                if (isLT)
                    if (dist ~= 0 && dist <= spatialS)
                        seed = seed * height;
                    end
                    if (t ~= 0)
                        seed = seed * height;
                    end
                end
                if (dist > spatialS)
                    seed = 0;
                end
                effectiveS = effectiveS + seed;
            end
        end
    end
    %disp(spatialS);
    %disp(temporalS);
    %disp(effectiveS);
end