%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kevin Griffin and Nathan Wei                                            %
% Active Grid Project: Data Presentation Script                           %
% Applies standard formatting to figures for use in paper preparation     %
% Operation: Run once with save_fig = 0, make edits to text and labels,   %
%            save as a .fig file manually, then run again with            %
%            save_fig = 1 to save the revised figure.                     %
% Dependencies: none                                                      %
% Last Edited:  08/31/2017                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------------------------- USER INPUTS ---------------------------- %%

fig_name    = 'corrf_lt3.9lt3_lt6.5lt5.fig'; % figure to modify
eps_name    = 'Fig747.eps'; % name to save .eps version under
fig_dir     = '../Figures/'; % relative or absolute path, where .fig lives
save_fig    = 1; % export as .eps and save revised .fig if 1

%% ---------------------------- FORMATTING ----------------------------- %%

close all; % make sure nothing else is open before starting

f1 = open(fullfile(fig_dir, fig_name)); % find and open figure

% Set parameters
set(f1,'color','w'); % figure must be in black and white
set(findall(f1,'-property','fontsize'),'fontsize',20); % must be legible
set(findall(f1,'-property','markersize'),'markersize',20);
set(findall(f1,'-property','linewidth'),'linewidth',2); % > 0.1 mm wide
set(findall(f1,'-property','interpreter'),'interpreter','latex');
f1_size = get(f1,'position');
f1_size(3:4) = [800 800]; % set new dimensions of figure
set(f1,'position',f1_size);

if save_fig
    cd(fig_dir);
    saveas(gcf, fig_name); % save revised figure under same name
    print(eps_name, '-depsc', '-tiff'); % export as eps
end