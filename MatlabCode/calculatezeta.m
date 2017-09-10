% Computes grid correlation shape parameter from a given 3D kernel.
% zeta is currently defined as the average difference between slopes of a
% perfect top hat and the given kernel (from outside in).
% Kevin Griffin and Nathan Wei, September 2017 (modified from CorrCheck.m)

function[zeta] = calculatezeta(sigmaS,sigmaT,height,shush)

% Don't print if told to shaddupyomouth
if nargin == 3 || ~shush
    if height == 0.1
        fprintf('lt%3.1flt%i  ', sigmaS, sigmaT)
    elseif height == 1
        fprintf('th%3.1fth%i  ', sigmaS, sigmaT)
    else
        fprintf('??%3.1f??%i  ', sigmaS, sigmaT)
    end
end

% Experiment parameters
meanU = 1.5; % m/s
timeStep = 0.1; % s
paddled = 0.115; % m

bound = ceil(sigmaS);
kernel = zeros(2*bound+1, 2*bound+1, 2*ceil(sigmaT)+1);
k_ref = zeros(2*bound+1, 2*bound+1, 2*ceil(sigmaT)+1);
% Build kernel
for j = -bound : bound
    for k = -bound : bound
        for t = -ceil(sigmaT) : ceil(sigmaT);
            dist = sqrt((j*j)+(k*k));
            if (dist <= 0)
                corrS = 1;
            elseif (dist <= sigmaS)
                corrS = height;
            else
                corrS = 0;
            end
            abs_t = abs(t);
            if (abs_t <= 0)
                corrT = 1;
            elseif (abs_t <= sigmaT)
                corrT = height;
            else
                corrT = 0;
            end
            kernel(j+bound+1,k+bound+1,t+ceil(sigmaT)+1) = corrS * corrT;
            if dist <= sigmaS && abs_t <= sigmaT
                k_ref(j+bound+1,k+bound+1,t+ceil(sigmaT)+1) = 1;
            end
        end
    end
end

% Calculate resulting paddle correlations
C_ker = convn(kernel, kernel);
C_ker = C_ker / max(C_ker(:)); % normalize by maximum value so max = 1
size_a = size(C_ker);
C_ref = convn(k_ref, k_ref);
C_ref = C_ref / max(C_ref(:)); % normalize by maximum value so max = 1

% measure of triangularity
zeta = 0;
vol = 0;
dV = paddled^2*timeStep*meanU;
cutoff = 1e-12;
h1 = (size_a(1)-1)/2+1;
h2 = (size_a(2)-1)/2+1;
h3 = (size_a(3)-1)/2+1;
corrLength = sum(C_ker(:)*dV/max(C_ker(:)))^(1/3);
for i = 1:size_a(1)
    for j = 1:size_a(2)
        for k = 1:size_a(3)
            % calculate distance to center element
            mydist = sqrt((i-h1)^2+(j-h2)^2+(meanU*timeStep/paddled)^2*(k-h3)^2);
            % compute difference btwn. perfect corr. slope and actual slope
            if C_ker(i,j,k) >= cutoff && corrLength-mydist ~= 0
                vol = vol + dV;
                % zeta = zeta + dV*(C_ref(i,j,k)-C_ker(i,j,k))/(mydist*((1/(1-C_ref(i,j,k)))-1)); % slope differences
                % zeta = zeta + dV*(C_ker(i,j,k))/(mydist*((1/(1-C_ref(i,j,k)))-1)); % just slope
                zeta = zeta + dV*(C_ker(i,j,k))/abs(mydist-corrLength); % l_alpha instead of rng - not bad w/ or w/o abs
                % GO BACK AND TRY USING L_ALPHA INSTEAD OF RNG BUT KEEP ALL
                % SLOPES POSITIVE. TRY THIS AND TRY WITH ACTUAL SLOPE AT
                % L_ALPHA (more physical but will probably fail)
            end
        end
    end
end

zeta = zeta / vol;

end