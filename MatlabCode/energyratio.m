dist = [9.6387 8.3454 7.2256 6.2561 5.4167 4.6899 4.0606 3.5158 3.0440 2.6356];
paddled = 0.115;

%weightedAverageEngergyArray for each kernel (lt5.2lt50 th3.9th3 lt5.2lt4)
E1 = [0.0027 0.0032 0.0037 0.0045 0.0050 0.0061 0.0075 0.0083 0.0103 0.0126];
E2 = [0.0021 0.0026 0.0030 0.0036 0.0042 0.0051 0.0056 0.0071 0.0083 0.0114];
E3 = [0.0021 0.0026 0.0030 0.0037 0.0044 0.0052 0.0064 0.0081 0.0099 0.0131];

%weightedAverageMASvsmArray for each kernel (lt5.2lt50 th3.9th3 lt5.2lt4)
U1 = [1.4423 1.4533 1.4512 1.4358 1.4277 1.4162 1.4035 1.3988 1.3721 1.3660];
U2 = [1.3943 1.3873 1.3781 1.3621 1.3629 1.3541 1.3251 1.3269 1.3156 1.3456];
U3 = [1.3918 1.3862 1.3828 1.3708 1.3677 1.3638 1.3667 1.3625 1.3544 1.3520];

e1 = E1./U1.^2;
e2 = E2./U2.^2;
e3 = E3./U3.^2;

b1= [0.1023, 8.2042, -1.0102];
b2= [0.0789, 9.6088, -0.9817];
b3= [0.0896, 10.9624, -1.0199];

x1 = linspace(10,100);
y1 = b1(1).*(x1-b1(2)).^b1(3);
y2 = b2(1).*(x1-b2(2)).^b2(3);
y3 = b3(1).*(x1-b3(2)).^b3(3);

figure(3);
hold on;
plot(x1,y1,'g');
scatter(dist/paddled,e1,1000,'.','g');
plot(x1,y2,'r');
scatter(dist/paddled,e2,1000,'.','r');
plot(x1,y3,'b');
scatter(dist/paddled,e3,1000,'.','b');
grid on;
h = gca;
set(h,'XScale','log');
set(h,'YScale','log');
set(h,'Fontsize', 20);
xlim([22.5 82]);
ylim([10^-3 7*10^-3]);
title('Power Law Fits');
ylabel('Normalized Energy (u^2 / U^2)');
xlabel('Normalized Distance (x / paddle width)');
legend('lt5.2lt50 160 minutes','lt5.2lt50 160 minutes','th3.9th3 80 minutes','th3.9th3 80 minutes','lt5.2lt4 80 minutes','lt5.2lt4 80 minutes');

figure(4);
hold on;
xlim([22.5 82]);
h = gca;
set(h,'Fontsize', 20);
set(h,'XScale','log');
plot(x1,y2./y1,'g');
scatter(dist/paddled,e2./e1,1000,'.','g');
plot(x1,y3./y1,'r');
scatter(dist/paddled,e3./e1,1000,'.','r');
plot(x1,y3./y2,'b');
scatter(dist/paddled,e3./e2,1000,'.','b');
title('Ratio of normalized energy');
legend('th3.9th3 / lt5.2lt50','th3.9th3 / lt5.2lt50','lt5.2lt4 / lt5.2lt50','lt5.2lt4 / lt5.2lt50','lt5.2lt4 / th3.9th3','lt5.2lt4 / th3.9th3');