%{ 

find power law fit. given a correlation function (smoothed)
f(r) = a*r(-m)
    

%}
L = hwils(MASC,sepval,2); %this is the integral length scale
%computes the array index start from which to fit 
for i=1 : length(sepval)
    if (sepval(i) > L) 
        start = i;
        break;
    end;
end;

figure(1);
hold on; 
[slope, intercept, MSE, R2] = logfit(sepvalc(start:length(sepvalc)),MASCc(start:length(MASCs))); 
%yApprox = (10^intercept)*sepval.^(slope);
hax = gca; 
line([L,L],get(hax,'YLim'), 'Color' , [1 0 0], 'LineWidth', 1); 
ylabel('correlation   ');
xlabel('distance (m)  ');

%this is the original function 

figure(2);
hold on; 
hax = gca; 
plot(sepvalc,MASCc,'LineWidth', 2)
line([L,L],get(hax,'YLim'), 'Color' , [1 0 0], 'LineWidth', 2); 
ylabel('correlation   ');
xlabel('distance (m)  ');




