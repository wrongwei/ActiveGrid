close all;
figure(1)
autocorr_kg(3.9, 3, true, 'c', 'm')
autocorr_kg(3.9, 3, false, 'b', 'r')
legend('th3.9th3','th3.9th3 corr','lt3.9lt3','lt3.9lt3 corr')

figure(2)
autocorr_kg(3.9, 3, true, 0, 'm')
autocorr_kg(6.5, 5, false, 0, 'b')
autocorr_kg(3.9, 3, false, 0, 'r')
legend('th3.9th3 corr', 'lt6.5lt5 corr', 'lt3.9lt3 corr')
%autocorr(5.2, 50, true, 0, 'c')
%autocorr(2, 50, false, 0, 'k')
%autocorr(5.2, 25, false, 0, 'g')