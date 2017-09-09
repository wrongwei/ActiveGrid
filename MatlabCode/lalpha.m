% Computes grid correlation parameters from a given 3D kernel.
% Kevin Griffin and Nathan Wei, September 2017 (modified from CorrCheck.m)

% MODIFY THIS -------------------------------------------------------------
function[corrLength] = lalpha(sigmaS,sigmaT,height)
% Kernel parameters
%sigmaS = 5.2;
%sigmaT = 25;
%height = 0.1;
if height == 0.1
    fprintf('lt%3.1flt%i  ', sigmaS, sigmaT)
elseif height == 1
    fprintf('th%3.1fth%i  ', sigmaS, sigmaT)
else
    fprintf('??%3.1f??%i  ', sigmaS, sigmaT)
end

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
        end
    end
end

% tmp, I am editing the kernel to make perfect kernels
% kernel_tmp = height*ones(size(kernel));
% kernel_tmp(bound+1,bound+1,ceil(sigmaT)+1) = 1;
% kernel = kernel_tmp;


% Calculate resulting paddle correlations
a_corr = convn(kernel, kernel);
a_corr = a_corr / max(a_corr(:)); % normalize by maximum value so max = 1
a_sum = sum(a_corr(:)); % integral of correlation
a_vol = (paddled^2*meanU*timeStep) * nnz(a_corr); % m^3
[aCs,~,aCt] = size(a_corr); % correlation lengths
a_avgcorr = max(a_corr(:)) / (a_vol^(1/3)); % normalize by volume of kernel
true_avg = a_sum/nnz(a_corr);
size_a = size(a_corr);
true_avg2 = a_sum/(size_a(1)*size_a(2)*size_a(3));

% measure of triangularity
triang = 0;
% rng = max dist of nonzero elm to center
rng = 0;
%cutoff = 0.001;
cutoff = 1e-12;
h1 = (size_a(1)-1)/2+1;
h2 = (size_a(2)-1)/2+1;
h3 = (size_a(3)-1)/2+1;
triCnt = 0;
for i = 1:size_a(1)
    for j = 1:size_a(2)
        for k = 1:size_a(3)
            if a_corr(i,j,k) >= cutoff
                mydist = sqrt((i-h1)^2+(j-h2)^2+(meanU*timeStep/paddled)^2*(k-h3)^2);
                if mydist > rng
                    rng = mydist;
                end
            end
        end
    end
end
corrLength = sum(a_corr(:)*paddled^2*timeStep*meanU/max(a_corr(:)))^(1/3);
for i = 1:size_a(1)
    for j = 1:size_a(2)
        for k = 1:size_a(3)
            % calc distance to center element
            mydist = sqrt((i-h1)^2+(j-h2)^2+(meanU*timeStep/paddled)^2*(k-h3)^2);
%             if a_corr(i,j,k) > cutoff
%                 triCnt = triCnt + 1;
%                 triang = triang + (rng - mydist) / a_corr(i,j,k);    
%             end

% my best guess
            if (rng - mydist) > 0 && a_corr(i,j,k) >= cutoff
                triCnt = triCnt + 1;
                triang = triang + a_corr(i,j,k) / (rng - mydist);
            end
            
%             if (2*corrLength - mydist) > 0 && a_corr(i,j,k) >= cutoff
%                 triCnt = triCnt + 1;
%                 triang = triang + a_corr(i,j,k) / (2*corrLength - mydist);
%             end
            
%             if (2*corrLength - mydist) > 0 && a_corr(i,j,k) >= cutoff
%                 triCnt = triCnt + 1;
%                 triang = triang + a_corr(i,j,k) / (2*corrLength - mydist);
%             end


%             if rng >= mydist & a_corr(i,j,k) >= cutoff
%                 triCnt = triCnt + 1;
%                 triang = triang - a_corr(i,j,k) + (rng - mydist)/rng;
%             end
%             if mydist > 0
%                 triCnt = triCnt + 1;
%                 triang = triang + (rng / mydist) * a_corr(i,j,k);
%             end
        end
    end
end
%true_avg3 = a_sum/n_cut;
%triangular_p = a_sum/nnz(a_corr)^2;

size(a_corr);
size_a_mod = [size_a(1),size_a(2),size_a(3)*meanU*timeStep/paddled];
%isotropy = (size_a(1)/size_a(2) + size_a(1)/(size_a(3)*meanU*timeStep/paddled) + size_a(2)/(size_a(3)*meanU*timeStep/paddled))/3
isotropy = min(size_a_mod)/max(size_a_mod);

%fprintf('sum_nnz %f  units %i nnz %i sum_ncut %f sum %f percent %f\n',a_sum/nnz(a_corr), size_a(1)*size_a(2)*size_a(3), nnz(a_corr), true_avg3, a_sum, nnz(a_corr)/(size_a(1)*size_a(2)*size_a(3)))
fprintf('(triang/triCnt) %7.2s\n',triang/triCnt)
%corrLength = (a_sum)^(1/3);

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
end