
% calibration file processing program
% HEADER

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