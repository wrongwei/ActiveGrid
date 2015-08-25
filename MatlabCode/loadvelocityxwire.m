%{
---------------------------------------------------------------------------
Data processing script for x-wire probes
Receives 2 data files of the same length (strings), a calibration file name
(string), a probe number (int), and a test temperature value (double) as
arguments
Returns vectors of streamwise (vx) and normal (vy) velocities, as well as
the number of points at which the calibration failed ('sorry')

Called from CheckIsoV2.m
Dependencies: processcalibrationxwire.m (which depends on an x-wire
calibration file and several velocity calibration .txt files),
extractprobeinfo.m, extractcalibdata.m, loadLabviewbin.m,
applytempcorrection.m

Nathan Wei, Willem van de Water, and Kevin Griffin
August 2015

NOTE: code currently operates on a 2nd-order polynomial fit for the
calibration curves. This is set in the calibration file, but is also
hard-coded into this program (rather than using the slow polyval function).
We do not recommend changing this, since higher-order polyfits could make
the calibration algorithm break down due to poor resemblance outside the
range of calibration.
---------------------------------------------------------------------------
%}

function [vx,vy,sorry] = loadvelocityxwire(file1, file2, calibrationfilename,...
    probenumber, actualtemp)

tstart = tic;
fprintf('  loading raw data from files: \n    %s\n    %s  \n',...
    file1, file2);

% Read the raw data
E1 = loadLabviewbin(file1);
E2 = loadLabviewbin(file2);

% Preventative error handling (we assume length(E1) = length(E2))
if length(E1) ~= length(E2)
    error('X-wire data set length mismatch for files %s and %s \n', file1, file2);
end

probeinfo = extractprobeinfo(calibrationfilename);

% Temperature correction
if ~isempty(actualtemp)
    fprintf('  calculating temperature correction...  \n');
    calibdata = extractcalibdata(calibrationfilename);
    % get resistances from calibration file
    Rc1 = probeinfo.probe(probenumber).sensor(1).resistance;
    Rw1 = Rc1 * (1 + probeinfo.probe(probenumber).sensor(1).overheat);
    Rc2 = probeinfo.probe(probenumber).sensor(2).resistance;
    Rw2 = Rc2 * (1 + probeinfo.probe(probenumber).sensor(2).overheat);
    % mean calibration temperature over all angles and velocities
    calibtemp = mean2(calibdata.calibtemp);
    E1 = applytempcorrection(E1, actualtemp, calibtemp, Rw1, Rc1, probeinfo.alpha);
    E2 = applytempcorrection(E2, actualtemp, calibtemp, Rw2, Rc2, probeinfo.alpha);
end

% Process calibration data
fprintf('  processing calibration data...  \n');
params1 = processcalibrationxwire(calibrationfilename, probenumber, 1);
params2 = processcalibrationxwire(calibrationfilename, probenumber, 2);
% params = n x p matrix (n is # of calib angles, p is order of polyfit)

% Velocity component vectors to return
vx = zeros(length(E1),1);
vy = zeros(length(E1),1);

% Set up parameters needed for calibration execution
vel1 = zeros(length(calibdata.angles),1);
vel2 = zeros(length(calibdata.angles),1);
diff = zeros(length(calibdata.angles),1);
sorry = 0;

%{
% Plot calibration curves for velocities vs. voltages (debugging)
range = linspace(0,10,100);
figure;
for j = 1 : length(calibdata.angles)
    plot(range,polyval(params1(j,:),range));
    hold on;
end
legend('1','2','3','4','5','6','7');
figure;
for j = 1 : length(calibdata.angles)
    plot(range,polyval(params2(j,:),range));
    hold on;
end
legend('1','2','3','4','5','6','7');
%}

% ------------------------- APPLY THE CALIBRATION -------------------------
fprintf('  converting voltages to velocities...  \n');
for i = 1 : length(E1)
    
    % use calibration curve fits and the given voltage to get one possible
    % velocity for every calibration angle
    for z = 1 : length(calibdata.angles)
        %vel1(z) = polyval(params1(z,:),E1(i)); % processcalibrationxwire polyfits
        %vel2(z) = polyval(params2(z,:),E2(i));
        % Using basic algebra instead of polyval is a lot faster for 2nd
        % order fits. Rewrite this if you use a higher-order fit (but why
        % would you want to do that anyway?)
        vel1(z) = params1(z,1)*(E1(i)^2)+params1(z,2)*(E1(i))+params1(z,3);
        vel2(z) = params2(z,1)*(E2(i)^2)+params2(z,2)*(E2(i))+params2(z,3);
        diff(z) = vel2(z) - vel1(z); % subtract two curves at each point
    end
    
    % find point of intersection (if there is one)
    t = 2;
    while t <= length(calibdata.angles)
        if (diff(t)*diff(t-1) <= 0) % if diff changes signs (or is 0)
            break; % t now represents point to the right of the intersection
        elseif (t == length(calibdata.angles))
            break; % reached end without finding anything :(
        else
            t = t + 1; % move on to the next angle if no intersection found
        end
    end
    % t now represents the right side of a possible range for the
    % intersection point (if no intersection found, t is at largest angle)
    % carry out linear interpolation btwn this point and its left neighbor
    angle = calibdata.angles(t-1) + ((vel2(t-1)-vel1(t-1))*...
        (calibdata.angles(t)-calibdata.angles(t-1))...
        / (vel1(t)-vel1(t-1)-vel2(t)+vel2(t-1)));
    velocity = ((vel2(t)-vel2(t-1))/(calibdata.angles(t)-...
        calibdata.angles(t-1)))*(angle-calibdata.angles(t-1))+vel2(t-1);
    %fprintf('Angle = %f, Velocity = %f \n',angle,velocity); % debugging
    
    % if the resulting angle is too big, try interpolating from first 2 points
    if (angle > max(calibdata.angles)*1.5)
        angle = calibdata.angles(2) + ((vel2(1)-vel1(1))*...
            (calibdata.angles(2)-calibdata.angles(1))...
            / (vel1(2)-vel1(1)-vel2(2)+vel2(1)));
        % if this angle is too small, then the calibration just doesn't
        % work here (curves don't intersect on a feasible range)
        if (angle < min(calibdata.angles)*1.5)
            sorry = sorry + 1; % number of points where calibration fails
            angle = NaN; % mark this point for future handling
        end
    end
    % produce components from calculated velocity and angle
    vx(i) = cos(angle * pi/180) * velocity;
    vy(i) = sin(angle * pi/180) * velocity;
end

totaltime = toc(tstart);
fprintf('  total time to load and convert the data: %.1f seconds.  \n', totaltime);
