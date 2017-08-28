% Plots a 2D representation of the correlation of a given set of angles, by
% radially averaging correlations about the given paddle.
% Nathan Wei, 12 August 2017 (modified from CorrCheck.m)

clear all;
close all;

% MODIFY THIS ------------------------------------------------------------------
A = load('/Users/weit/Documents/nwei/College/Engineering/MAE 339/angles_th3.9th3_0813_1.txt'); %load in the angle data file generated from menuII (with the angle header removed)
B = load('/Users/weit/Documents/nwei/College/Engineering/MAE 339/angles_lt3.9lt3_h0.1_0811.txt');
%choose paddle with which you want to see the correlation of other paddles 
padrow = 7;
padcol = 7;

% ---------------------------------------------------------------------
sz_A = size(A);
time = sz_A(1);
norm = 0;
a = zeros(11,13); % "a" is the row*col matrix. a(i,j) is the correlation of that paddle with the (padrow, padcol) paddle.
dist = zeros(11,13); % distance of each paddle from chosen paddle

%calculate correlation for each time slice
for t = 1 : time
    timeslice = A(t,:);
    agrid = reshape(timeslice,13,11);
    agrid = transpose(agrid); %convert timeslice to 11*13 matrix
    norm = norm + agrid(padrow,padcol)*agrid(padrow,padcol); %normalization
    
    for i = 1 : 11
        for j = 1 : 13
            a(i,j) = a(i,j) + agrid(i,j)*agrid(padrow,padcol);
        end
    end
end

a = a./time; % calculate average correlation across all time slices
a = a./(norm/time); % normalization

%calculate distance to every other paddle (no wrapping around)
for i = 1 : 11
    for j = 1 : 13
        dist(i,j) = sqrt((i-padrow)^2+(j-padcol)^2);
    end
end

% Calculate average value at each distance
dists = unique(dist);
a_rad = zeros(length(dists), 1);
for i = 1 : length(dists)
    ind = find(ismember(dist, dists(i))); % find equidistant points (linear indices)
    a_rad(i) = mean(a(ind)); % average all points at given distance
end

% ----- DO IT AGAIN FOR THE SECOND FILE -----
sz_B = size(B);
timeb = sz_B(1);
normb = 0;
b = zeros(11,13); % "a" is the row*col matrix. a(i,j) is the correlation of that paddle with the (padrow, padcol) paddle.
distb = zeros(11,13); % distance of each paddle from chosen paddle

%calculate correlation for each time slice
for t = 1 : timeb
    timeslice = B(t,:);
    bgrid = reshape(timeslice,13,11);
    bgrid = transpose(bgrid); %convert timeslice to 11*13 matrix
    normb = normb + bgrid(padrow,padcol)*bgrid(padrow,padcol); %normalization
    
    for i = 1 : 11
        for j = 1 : 13
            b(i,j) = b(i,j) + bgrid(i,j)*bgrid(padrow,padcol);
        end
    end
end

b = b./timeb; % calculate average correlation across all time slices
b = b./(normb/timeb); % normalization

%calculate distance to every other paddle (no wrapping around)
for i = 1 : 11
    for j = 1 : 13
        distb(i,j) = sqrt((i-padrow)^2+(j-padcol)^2);
    end
end

% Calculate average value at each distance
distsb = unique(distb);
b_rad = zeros(length(distsb), 1);
for i = 1 : length(distsb)
    ind = find(ismember(distb, distsb(i))); % find equidistant points (linear indices)
    b_rad(i) = mean(b(ind)); % average all points at given distance
end

% % plot correlation vs. distance
figure();
d = dist(:);
corr = a(:);
scatter(d,corr,10,'bo','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on;
plot(dists, a_rad, 'bs', 'MarkerSize', 10);

db = distb(:);
corrb = b(:);
scatter(db,corrb,10,'ro','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on;
plot(distsb, b_rad, 'rs', 'MarkerSize', 10);

l1 = legend('Top hat', 'Top hat (radial-averaged)', ...
    'Long tail', 'Long tail(radial-averaged)');
set(l1, 'FontSize', 14);
set(l1, 'interpreter', 'latex');

xlabel('Radial Distance (paddles)', 'interpreter', 'latex', 'fontsize', 16);
ylabel('Time-Averaged Normalized Correlation', 'interpreter', 'latex', 'fontsize', 16);
title('Time-Averaged Correlation vs. Radial Distance for Two Kernels', ...
    'fontsize', 20, 'interpreter', 'latex');
ylim([0,1]);