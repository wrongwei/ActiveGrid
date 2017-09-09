%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kevin Griffin and Nathan Wei                                            %
% Active Grid Project: Data Analysis Script                               %
% Compute flow properties of the tunnel with grid open                    %
% Dependencies: none                                                      %
% Last Edited:  09/03/2017                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate mean velocity variance, Reynolds stresses (?), baseline
% turbulence intensity, etc.

%% ---------------------------- USER INPUTS ---------------------------- %%
path = '/Users/weit/Documents/nwei/College/Engineering/MAE 339/data_0828';
file = 'statscorr_opengrid_0828.mat';
rho  = 1.195; % kg/m^3

%% ---------------------------- PROCESSING ----------------------------- %%

% Get velocities from .mat file
data = load(fullfile(path, file), 'u', 'deltaT');
u = data.u(1:ceil(end)); % velocity signal, in m/s
dT = data.deltaT;
clear data;
u_f = u - mean(u); % fluctuation velocities

% Calculate mean velocity variance
u_var = var(u); % filter this?

% Calculate Reynolds stresses: TauR_i,j = -rho mean(ui'uj')
tau_Rxy = rho * mean(u_f.^2); % Pa
% Note that i and j are different components of the velocity. If we assume
% isotropic turbulence, ui' = uj' for all three flow dimensions.

% Calculate turbulence intensity (percent)
u_int = (rms(u_f) / mean(u)) * 100;

%% --------------------------- PLOT RESULTS ---------------------------- %%

% Print results
fprintf('Variance in mean velocity: %.4f m/s \n', u_var);
fprintf('Average Reynolds stress (x-y): %.4f Pa \n', tau_Rxy);
fprintf('Turbulence intensity: %.4f %% \n', u_int);

%{
% Plot velocity signal
t = (0 : 1 : length(u) - 1) * dT;
figure();
plot(t, u);
xlabel('Time (s)');
ylabel('Velocity (m/s)');

% Plot Reynolds stresses
figure();
plot(t, tau_Rxy);
xlabel('Time (s)');
ylabel('Reynolds Stress (Pa)');
%}