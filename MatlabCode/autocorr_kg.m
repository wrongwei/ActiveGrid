function [  ] = autocorr_kg( spatialSig, temporalSig, tophat, plotStrKer, plotStrCorr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%end
%clear all;
%close all;
% arguments:
% spatialSig is the spatial sigma
% temporalSig is the temporal sigma
% tophat is a boolean that is true for tophat kernels and false for longtails
% then there are two formating strings for the plots
% use 0 if no plot is desired
% use ' ', ' ' if no formatting is desired

tol = 1e-10;
if abs(1.3 - spatialSig/temporalSig) < tol
    fprintf('This kernel is isotropic! Proceeding...\n')
else
    fprintf('Anisotropic kernels not supported!!!\n')
    assert(0) %Anisotropic kernels not supported
end
s = temporalSig;
h = 0.1;
alpha = 1;
kernel = ones(1+2*s,1);
xKer = zeros(1+2*s,1);
for i = 1:1+2*s
    xKer(i) = -s-1+i;
end
if tophat
    fprintf('top hat detected\n');
else
    fprintf('long tail assumed\n');
    % modify the kernel to make it a longtail
    for i = 1:s-alpha
        kernel(i) = h;
        kernel(i + 1 + alpha + s) = h;
    end
end

%perform the convolution of kernel with itself
%corrL = length(kernel);
%longxaxis = xaxis;
width = 2*length(kernel)+1;
xVals = zeros(width,1);
for i = 1:width
    xVals(i) = 2*(-s-1)+i;
end

result = zeros(width,1);
for mainCnt = 1:width
    for stpCnt = mainCnt-s:mainCnt+s
        mainCnt
        xVals(mainCnt)
        stpCnt
        % check if in bounds of stamp kernel and within s units of main kernel
        if all([stpCnt >= 1, stpCnt <= width])
            if all([xVals(stpCnt) <= s, xVals(stpCnt) >= -s])
                % the result += contribution from stamp kernel * that of main kernel
                result(mainCnt) = result(mainCnt) + kernel(stpCnt+s+1-mainCnt)*kernel(stpCnt-s-1)
            end
        end
    end
end
% for stpCnt = (1+length(kernel)xVals(mainCnt)):(length(result)+xVals(mainCnt))
%     skip the part of the stamping kernel that is left of the main kernel
%     if stpCnt < 1
%         result(i) += 0;
%         skip the part of the stamping kernel that is right of the main kernel
%     elseif stpCnt > length(kernel)
%         result(i) += 0;
%     else
%         xVals(mainCnt)
%         result(mainCnt) = result(mainCnt) + kernel(stpCnt-xVals(mainCnt))*kernel(mainCnt)
%     end
% end
% end
result = result./max(result)
%close all;
%for plotting we want the kernels to have vertical lines
if ~tophat
    %for long tails additional points must be added to the kernel to get
    %vertical lines at the alpha core
    xKerPlot = zeros(length(xKer)+2,1);
    xKerPlot(1) = xKer(1); xKerPlot(end) = xKer(end);
    kernelPlot = zeros(length(kernel)+2,1);
    kernelPlot(1) = 0; kernelPlot(end) = 0;
    j = 1;
    for i = 1:length(xKer)
        kernelPlot(j) = kernel(i);
        xKerPlot(j) = xKer(i);
        if (i == s-1)
            j = j + 1;
            xKerPlot(j) = xKer(i+1);
            kernelPlot(j) = kernel(i);
        elseif (i == s + 2*alpha)
            j = j + 1;
            xKerPlot(j) = xKer(i);
            kernelPlot(j) = kernel(i+1);
        end
        j = j + 1;
    end
    xKer = xKerPlot
    kernel = kernelPlot
end
xKerPlot = zeros(length(xKer)+2,1);
xKerPlot(1) = xKer(1); xKerPlot(end) = xKer(end);
kernelPlot = zeros(length(kernel)+2,1);
kernelPlot(1) = 0; kernelPlot(end) = 0;
for i = 1:length(xKer)
    xKerPlot(i+1) = xKer(i);
    kernelPlot(i+1) = kernel(i);
end

if (plotStrKer ~= 0)
    assert(ischar(plotStrKer)) % should be enclosed in single quotes for example: 'd'
    plot(xKerPlot, kernelPlot, plotStrKer)
    hold on;
    legend('Kernel')
end
if (plotStrCorr ~= 0)
    assert(ischar(plotStrCorr)) % should be enclosed in single quotes for example: 'd'
    plot(xVals, result, plotStrCorr)
    hold on;
    legend('Autocorrelation')
end
if (plotStrKer ~= 0 & plotStrCorr ~=0)
    legend('Kernel','Autocorrelation')
end

end