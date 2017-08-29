% Plots a 2D representation of the correlation of a given set of angles, by
% radially averaging correlations about the given paddle.
% Nathan Wei, 28 August 2017 (modified from CorrCheck.m)

% MODIFY THIS -------------------------------------------------------------

% Kernel parameters
sigmaS = 6.5;
sigmaT = 5;
height = 1;

% Experiment parameters
meanU = 1.5; % m/s
timeStep = 0.1; % s
paddled = 0.115; % m

% -------------------------------------------------------------------------
bound = ceil(sigmaS);
kernel = zeros(2*bound+1, 2*bound+1, 2*ceil(sigmaT)+1);
% Build kernel
for j = -bound : bound
    for k = -bound : bound
        for t = -ceil(sigmaT) : ceil(sigmaT);
            dist = sqrt((j*j)+(k*k));
            if (dist == 0)
                corrS = 1;
            elseif (dist <= sigmaS)
                corrS = height;
            else
                corrS = 0;
            end
            abs_t = abs(t);
            if (abs_t == 0)
                corrT = 1;
            elseif (abs_t <= sigmaT)
                corrT = height;
            else
                corrT = 0;
            end
            kernel(j+bound+1,k+bound+1,t+ceil(sigmaT)+1) = corrS * corrT;
        end
    end
end
% Calculate resulting paddle correlations
a_corr = convn(kernel, kernel);
% a_corr = a_corr / max(a_corr(:)); % normalize by maximum value so max = 1
a_sum = sum(a_corr(:)); % integral of correlation, minus central paddle
a_vol = (paddled^2*meanU*timeStep) * nnz(a_corr); % m^3
[aCs,~,aCt] = size(a_corr); % correlation lengths
a_avgcorr = max(a_corr(:)) / (a_vol^(1/3)); % normalize by volume of kernel
% a_length = a_vol^(1/3) % characteristic correlation length for kernel
%{
% plot correlation vs. distance
figure();
hold on;
% plot(dists(:), a_corr(:), 'r.');
surf(a_corr(:,:,end/2+0.5));

xlabel('Radial Distance, $\frac{r}{\sigma_e}$', 'interpreter', 'latex', 'fontsize', 16);
ylabel('Time-Averaged Normalized Correlation', 'interpreter', 'latex', 'fontsize', 16);
title('Time-Averaged Correlation vs. Radial Distance for Two Kernels', ...
    'fontsize', 20, 'interpreter', 'latex');
%ylim([0,1]);
%}