% plots the mesh, contour plot and 3D scatter plot of correlation of all
% other paddles on the grid with a user specified paddle (row, col)
% in the 3D plots the color is the correlation 
% Used to test accuracy of kernels written in C++
% Horace Zhang + Jessie Liu Summer 2014
% Modified by Nathan Wei, Summer 2015

clear all;
close all;

% MODIFY THIS ------------------------------------------------------------------
pathstring = fileparts('/Users/kevin/Documents/Data/data08_07_15/');
%pathstring = fileparts('/Users/nathan/Documents/Data/data07_30_15/');
addpath(pathstring);
A = load('angles_lt1.3lt1_h0.4_Re_0807.txt'); %load in the angle data file generated from menuII (with the angle header removed)
%choose paddle with which you want to see the correlation of other paddles 
padrow = 7;
padcol = 7;

% ---------------------------------------------------------------------
size = size(A);
time = size(1);
norm = 0;
a = zeros(11,13); % "a" is the row*col matrix. a(i,j) is the correlation of that paddle with the (padrow, padcol) paddle.
padangle = zeros(time,1); % stores angle at each timeslice for chosen paddle
inst_corr = zeros(11,13,time); % stores correlation of each paddle with chosen paddle, over time
avgangle = zeros(time,1); % stores average angle of each timeslice
rmsd = zeros(time,1); % stores rmsd of each timeslice

%calc avg angle value
avgAllAngles = 0;
for t = 1 : time
    timeslice = A(t,:);
    agrid = reshape(timeslice,13,11);
    agrid = transpose(agrid); %convert timeslice to 11*13 matrix
    for i = 1 : 11
        for j = 1 : 13
            avgAllAngles = avgAllAngles + abs(agrid(i,j));
        end
    end
end
avgAllAngles = avgAllAngles / time / 13.0 / 11.0

%calc avg rms
rmsAllAngles = 0;
for t = 1 : time
    timeslice = A(t,:);
    agrid = reshape(timeslice,13,11);
    agrid = transpose(agrid); %convert timeslice to 11*13 matrix
    for i = 1 : 11
        for j = 1 : 13
            rmsAllAngles = rmsAllAngles + (abs(agrid(i,j))-avgAllAngles)^2;
        end
    end
end
rmsAllAngles = sqrt(rmsAllAngles / time / 13.0 / 11.0)

for t = 1 : time
    timeslice = A(t,:);
    agrid = reshape(timeslice,13,11);
    agrid = transpose(agrid); %convert timeslice to 11*13 matrix
    norm = norm + (abs(agrid(padrow,padcol))-avgAllAngles)*(abs(agrid(padrow,padcol))-avgAllAngles); %normalization
    padangle(t) = agrid(padrow,padcol);
    
%calculate correlation
    for i = 1 : 11
        for j = 1 : 13
            a(i,j) = a(i,j) + (abs(agrid(i,j))-avgAllAngles)*(abs(agrid(padrow,padcol))-avgAllAngles);
            %inst_corr(i,j,t) = (agrid(i,j)*agrid(padrow,padcol))/(agrid(padrow,padcol)*agrid(padrow,padcol));
            avgangle(t) = avgangle(t) + (agrid(i,j) / 143);
        end
    end
end

a = a./time; % calculate average
a = a./(norm/time); %normalization
%{
%calculate temporal correlation (CURRENTLY NONFUNCTIONAL)
temp_norm = 0;
padcorr = zeros(time,1); % stores unnormalized correlation for each timeslice

for i = 1 : time
    temp_norm = temp_norm + (padangle(i)*padangle(i)); % same job as norm
    for j = 1 : time
        if i == j
            continue;
        else
            padcorr(i) = padcorr(i) + (padangle(i)*padangle(j));
        end
    end
end

padcorr = padcorr./(time - 1); % calculate average over time-1 comparisons
padcorr = padcorr./(temp_norm/(time - 1)); % normalization
%}
%calculate distance to every other paddle (no wrapping around)

for i = 1 : 11
    for j = 1 : 13
        dist(i,j) = sqrt((i-padrow)^2+(j-padcol)^2);
    end
end

for t = 1 : time
    agrid = transpose(reshape(A(t,:),13,11));
    for i = 1 : 11
        for j = 1 : 13
            rmsd(t) = rmsd(t) + (agrid(i,j)-avgangle(t))^2 / 143;
        end
    end
    rmsd(t) = sqrt(rmsd(t));
end

%plot mesh
figure(1);
set(gca, 'fontsize', 45)
m2 = meshc(a);
set(m2, 'LineWidth', 2)
xlabel('Columns');
ylabel('Rows');
colorbar;
caxis([0 1]); %keep the colorbar scale the same

%plot contour map
figure(2);
set(gca, 'fontsize', 45)
m3 = contour(a, 'LineWidth', 2);
grid on;
xlabel('Columns');
ylabel('Rows');
set(gca,'XTick',1:13);
set(gca,'YTick',1:11);
caxis([0 1]);
colorbar;

% plot correlation vs. distance
figure(3);
d = dist(:);
corr = a(:);
scatter(d,corr);
xlabel('Distance (paddles)');
ylabel('Correlation');
ylim([0,1]);

fid = fopen('res.txt','wt');
for i=1:length(d)
     fprintf(fid,'%.4f  %.4f  \n', d(i), corr(i));
end;

% plot temporal correlation of chosen paddle vs. time
%{
figure(4);
plot(padcorr,'o');
xlabel('Time units'); % temporal spacing between timeslices
ystring = sprintf('Correlation of Paddle (%d,%d)',padcol,padrow);
ylabel(ystring);
%}
% plot RMSD of grid over time
figure(4);
plot(rmsd);
title('RMSD of grid vs. time');
xlabel('Time units');
ylabel('Angle (degrees)');
fprintf('Average RMS = %.4f \n', mean(rmsd));

%{
% Angle visualizer (shows changes over time and saves result as movie)
disp(time);
for i = 1 : time;
    figure(5);
    disp(i);
    set(gca, 'fontsize', 45)
    m5 = meshc(transpose(reshape(A(i,:),13,11)));
    set(m5, 'LineWidth', .3)
    title('Paddle angles over time');
    xlabel('Columns');
    ylabel('Rows');
    zlim([-90 90]);
    caxis([-90 90]);
    colorbar;
    F(i) = getframe;
    drawnow
end
% Play movie (frames array, number of times, frame rate)
%movie(F, 2)

% Save movie
writerObj = VideoWriter('test', 'MPEG-4'); % create VideoWriter object
writerObj.FrameRate = 30; % set frame rate (default = 30)
open(writerObj); % open VideoWriter object for editing
for i = 1 : length(F);
    writeVideo(writerObj, F(i)); % write in frames one by one
end
close(writerObj); % close VideoWriter object (finishes movie-making)
disp('Saved test.mp4!');
%}

