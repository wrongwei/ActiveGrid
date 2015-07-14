% Uses nlinfit to determine coefficients in power law fit 
% y = a(x-b)^c
% U is the mean velocity vector 
% eps is the energy vector u^2
% d is the distances vector (meters) 
% b0 = [a0,b0,c0] intial paramters for nonlinear power law fit
% outputs are the coefficients, the R^2  value of the fit, and a 95% 
% confidence interval on the value of c (the exponent in the power law)
% plots the fit on normal and loglog axes
% Horace Zhang + Jessie Liu Summer 2014
% Dependencies: none 
function [a, b, c, R2, conf] = Edecfit(dist, U, eps, b0)
    % xscale for plotting
    paddled = 0.13; % distance between two adjacent paddles (m)
    %normalize
    normeps = eps./U.^2;
    normd = dist/paddled;
    xscale = min(normd): 0.1: max(normd);
    
        % uses nlinfit to be original data
        func = @(b,x)b(1).*(x-b(2)).^b(3); %the power law y = a(x-b)^c
        %initb = [1, 3, -1]; 
        initb = b0; 
        [b1, R, J, covb, mse] = nlinfit(normd,normeps, func, initb); 
        CI = nlparci(b1, R, 'covar', covb);

        % this is roughly the plus/minus variation away from the mean (half the interval): 
        exponent_error = diff(CI(3,:))/2; 
        
        %the fit on normal axes
        figure(1);
        hold on;
        grid on;
        h = gca;
        set(h,'XScale','linear');
        set(h,'YScale','linear');
        set(h,'Fontsize', 20);
        scatter(normd, normeps, 1000, 'b', '.')
        plot(xscale, b1(1).*(xscale - b1(2)).^b1(3), 'b' )
        ylabel('Normalized Energy   '); 
        xlabel('Normalized Distance '); 
        
        
        % the fit on loglog axes
        figure(2); 
        hold on;
        grid on;
        h = gca;
        set(h,'XScale','log');
        set(h,'YScale','log');
        set(h,'Fontsize', 20);
        scatter(normd, normeps, 1000, 'b', '.')
        plot(xscale, b1(1).*(xscale - b1(2)).^b1(3), 'b' )
        ylabel('Normalized Energy   '); 
        xlabel('Normalized Distance '); 
       
        
        a = b1(1); 
        b = b1(2);
        c = b1(3);
        SStot = sum((normeps - mean(normeps)).^2);
        SSres = sum(R.^2);
        R2 = 1 - SSres/SStot;
        % this is the 95% Confidence Interval on the exponent: c 
        % plus/minus this number is the 95% confidence interval
        conf =  exponent_error / -b1(3);
        fprintf('[a,b,c,R2,conf] is [%.4f, %.4f, %.4f, %.4f, %.4f] \n', [a,b,c,R2,conf]);
        fprintf('The 95%% confidence range is:  [%.4f, %.4f, %.4f] \n', [c-conf, c, c+conf]); 
   
   
end