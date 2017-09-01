%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nathan Wei and Kevin Griffin                                            %
% Active Grid Project: Data Presentation Script                           %
% Applies standard formatting to figures for use in paper preparation     %
% Dependencies: export_fig.m?                                             %
% Last Edited:   08/30/2017                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------------------- USER INPUTS ---------------------------- %%

fig_name    = 'corrf_lt3.9lt3_lt6.5lt5.fig';
fig_dir     = '../Figures/'; % relative or absolute path
save_fig    = 1; % export as .eps and resave .fig if 1

exfigpath   = '/Users/weit/Documents/MATLAB'; % location of export_fig.m

% NOTE: To export to .eps, you need to have GhostScript and Xpdf installed.
% See https://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% NOTE 2: Currently using print to export to .eps since I can't install the
% extensions for export_fig.

%% ---------------------------- PROCESSING ----------------------------- %%

f1 = open(fullfile(fig_dir, fig_name)); % find and open figure

% Set parameters
set(f1,'color','w');
set(findall(f1,'-property','fontsize'),'fontsize',20); % was 24
set(findall(f1,'-property','interpreter'),'interpreter','latex');
f1_size = get(f1,'position');
f1_size(3:4) = [1000 1000]; % set new dimensions of figure
set(f1,'position',f1_size);

if save_fig
    cd(fig_dir);
    addpath(exfigpath); % access to export_fig ensemble
    saveas(gcf, fig_name);
    % export_fig(fig_name(1:(end-4)), '-nocrop', '-r300', '-pdf');
    fig_name_new = [fig_name(1:(end-4)), '.eps'];
    print(fig_name_new, '-depsc', '-tiff'); % export as eps
    rmpath(exfigpath);
end