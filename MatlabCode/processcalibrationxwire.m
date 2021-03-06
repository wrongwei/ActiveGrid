%{
---------------------------------------------------------------------------
Calibration processing script for x-wire probes
Receives a calibration file name (string), a probe number (int), and a
sensor number (int) as arguments
Returns a N x P matrix with polynomial fit coefficients
 - N is the number of calibration angles (number of velocity calibration
    data sets, as referenced from the x-wire calibration file)
 - P is the order of the polynomial fit + 1 (currently using 2nd order fit,
    so P = 3)

Called from loadvelocityxwire.m
Dependencies: x-wire calibration file and N velocity calibration .txt files
(generated by StreamWare, w/ 9 columns) should be present in the directory

Nathan Wei, Willem van de Water, and Kevin Griffin
August 2015
---------------------------------------------------------------------------
%}

function coefficients = processcalibrationxwire(calibrationfilename, ...
    probenumber, sensornumber)

% ---- load probeinfo, calibdata, ...
[pathstr, name] = fileparts(calibrationfilename);
odiry = pwd;
if (~isempty(pathstr)), cd(pathstr); end
run(name);
cd(odiry);

% ---- organize the data:
% maximum channel number:
maxchannel = 0;
for i = 1:probeinfo.numprobes
    for j = 1:probeinfo.probe(i).numsensors
        if (probeinfo.probe(i).sensor(j).channel > maxchannel)
            maxchannel = probeinfo.probe(i).sensor(j).channel;
        end
    end
end

coefficients = zeros(length(calibdata.angles), calibdata.polyfit + 1);
for n = 1 : length(calibdata.angles)
    % not sure what this does, but leaving it here commented out
    %{
    if (size(calibdata.v, 1) == 1)
        calibdata.v = repmat(calibdata.v, [maxchannel 1]);
    end
    if (size(calibdata.T, 1) == 1)
        calibdata.T = repmat(calibdata.T, [maxchannel 1]);
    end
    %}
    
    % ---- sort data by increasing velocity:
    %{
    channel = probeinfo.probe(probenumber).sensor(sensornumber).channel;
    [calibdata.v(n, channel, :), I] = sort(calibdata.v(n, channel, :));
    if ~isempty(calibdata.fanspeed),   calibdata.fanspeed = calibdata.fanspeed(I);   end
    if ~isempty(calibdata.dp),   calibdata.dp = calibdata.dp(I);   end
    if ~isempty(calibdata.dens) && length(calibdata.dens) > 1
        calibdata.dens = calibdata.dens(I);
    end
    if (length(calibdata.T(channel, :)) > 1)
        calibdata.T(channel, :) = calibdata.T(channel, I);
    end
    calibdata.E(n, channel, :) = calibdata.E(n, channel, I);
    %}
    % the relevant data:
    thisV = calibdata.v(n, sensornumber, :);
    thisE = calibdata.E(n, sensornumber, :);
    % polyfit returns a vector of coefficients for fit of order 'polyfit'
    coefficients(n, :) = polyfit(thisE, thisV, calibdata.polyfit);
end