
% distill the data: dissipation, spectrum, ESS, ...  
% for single-wire data.  
% 
% uses: 
% rs, deltaT, vsm, vss, 
% S2, S3, S4, S6, 
% eps, 
% freq, specf

loeta = 20; 
upeta = 200; 

singleprobe = false; 
  % if singleprobe is false, give probe name: 
probename = 'S1'; 

odiry = pwd; 
warning off

% ---- work on data in directory: 
%{
diry = '/Users/gregory/Desktop/wind tunnel data/1 bar SF6/Sep 5 NSTAP 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/1 bar SF6/Sep 5 NSTAP 20 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/2 bar SF6/Nov 11 NSTAP 12.5 Hz'; 
diry = '/Users/gregory/Desktop/wind tunnel data/4 bar SF6/Aug 14 NSTAP III 20 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/8 bar SF6/Aug 20 NSTAP II 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/12 bar SF6/Aug 21 NSTAP 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/14 bar SF6/Oct 14 NSTAP 12.5 Hz'; 
diry = '/Users/gregory/Desktop/wind tunnel data/15 bar SF6/Aug 22 NSTAP 12.5 Hz/'; 

diry = '/Users/gregory/Desktop/wind tunnel data/1 bar Air/Nov 24 structure 12.5 Hz/'; 
%diry = '/Users/gregory/Desktop/wind tunnel data/1 bar SF6/Sep 2 profile and structure 20 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/1 bar SF6/Sep 16 structure 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/2 bar SF6/Sep 19 structure 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/4 bar SF6/Sep 20 structure 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/8 bar SF6/Nov 2 structure 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/12 bar SF6/Oct 21 structure 12.5 Hz/'; 
diry = '/Users/gregory/Desktop/wind tunnel data/14 bar SF6/Oct 18 structure 12.5 Hz'; 
%diry = diries{i}; 
%}
diry = '/Users/Horace/documents/Germany2014/MATLABCode/morecode';
cd(diry)
load reduceddata
%load reducedata
run('params'); 
fluidinfo = extractfluidinfo(calibname); 
cd(odiry)
nu = fluidinfo.nu; 

if ~singleprobe
	renamevars; 
end

  % basic stats: 
fprintf(diry); 
fprintf('\n'); 
fprintf('mean flow speed:        %.2f m/s \n', mean(vsm(:))); 
fprintf('turbulence intensity:   %.2f %% \n', 100*mean(vss(:))/mean(vsm(:))); 

  % dissipation from 3rd order structure: 
ra = rs*deltaT*mean(vsm(:)); 
eps3 = max((S3./ra')*(5/4)); 
fprintf('dissipation from S3:    %.4f m^2/s^3 \n', eps3); 

  % normalized dissipation: 
L = hwils(1-S2./(2*mean(vss(:))^2), ra); 
meandiss = mean(mean(eps(:,:,1)));  % where does this come from?  
meandiss = mean(eps(:,1)); 
normdiss = meandiss*L/mean(vss(:))^3; 
fprintf('integral scale:         %.1f mm \n', 1000*L); 
fprintf('dissipation rate:       %.4f m^2/s^3 \n', meandiss); 
fprintf('epsilon L / u^3:        %.2f \n', normdiss); 

%figure; 
%loglog(ra, S3./(meandiss*ra'), 'kx')
%hold on
%grid, grid minor
%plot([min(ra) max(ra)], 0.8*[1 1], 'k-')
%xlabel('r [m]')
%ylabel('S_3 / \epsilonr')

  % Reynolds number: 
lambda = sqrt(15*nu/meandiss)*mean(vss(:)); 
Rlambda = mean(vss(:))*lambda/nu; 
fprintf('Taylor scale:           %.2f mm \n', 1000*lambda); 
fprintf('Reynolds number:        %d \n', round(Rlambda)); 

  % make the spectrum: 
distillspectrum; 
fprintf('Kolmogorov scale:       %d microns \n', round(10^6*eta)); 

figure, loglog(Xspec/eta, Yspec)
hold on, grid, grid minor
title(diry)
xlabel('k_1 [1/m]')
ylabel('E_{11}(k_1) [m^2/s^2]')

figure, loglog(Xspec, Yspeck)
hold on, grid, grid minor
title(diry)
xlabel('k_1 \eta')
ylabel('E_{11}(k_1) / \epsilon^{2/3}k_1^{-5/3}')

  % make ESS plots: 
  % find scaling exponent between 10 and 100 eta: 
ll = max(find(ra/eta < loeta)); 
ul = min(find(ra/eta > upeta)); 
zeta2 = mean( diff(log(S2(ll:ul))) ./ diff(log(S3a(ll:ul))) ); 
zeta4 = mean( diff(log(S4(ll:ul))) ./ diff(log(S3a(ll:ul))) ); 
zeta6 = mean( diff(log(S6(ll:ul))) ./ diff(log(S3a(ll:ul))) ); 
fprintf('2nd order ESS exponent: %0.3f \n', zeta2); 
fprintf('4nd order ESS exponent: %0.3f \n', zeta4); 
fprintf('6nd order ESS exponent: %0.3f \n', zeta6); 

%figure, semilogx(ra(2:end)/eta, diff(log(S2))./diff(log(S3a)), 'kx')
%hold on, grid, grid minor
%semilogx(ra(2:end)/eta, diff(log(S4))./diff(log(S3a)), 'kx')
%semilogx(ra(2:end)/eta, diff(log(S6))./diff(log(S3a)), 'kx')
%axis([10^0 10^3 0 2.2])
%title(diry)
%xlabel('r / \eta')
%ylabel('d log(S_n) / d log(S_3)')


warning on
