%Correlation Clipper is a handy file for storing a correlation function, so
%that it takes up less memory and be sent by email to someone that doesn't
%have matlab

number_of_points = 200;
length_of_signal = length(sepvalc);
%throw out the remainder
spacing = floor(length_of_signal/number_of_points);

smallCorr3 = zeros(number_of_points,1);
smallSepvalc3 = zeros(number_of_points,1);

for i = 1:length_of_signal
    
   if (mod(i, spacing) == 0)
       smallSepvalc3(i/spacing) = sepvalc(i);
       smallCorr3(i/spacing) = MASCc(i);
   end
end
figure(1);
set(gca,'yscale','log');
plot(smallSepvalc3,smallCorr3);
hold on;
plot(sepvalc,MASCc);
