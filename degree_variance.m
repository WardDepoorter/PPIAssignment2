clear;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);


%% Load and prepare observed data (OBS)
nmax = 180;
Vmax = (nmax+1)*(nmax+2)/2;

OBS_full = zeros(Vmax,4);

OBS_raw = readmatrix('GRAILdata2.txt'); % n,m,C,S
OBS_data = OBS_raw(:,1:4); 

% First row is C0,0 = 0
OBS_full(1,:) = [0 0 0 0];

% Loop to fill remaining rows safely
% k counts the position in OBS_full
k = 2;
for idx = 1:size(OBS_data,1)
    n = OBS_data(idx,1);
    m = OBS_data(idx,2);
    if k > Vmax
        break; % don't exceed expected size
    end
    OBS_full(k,1:2) = [n m];
    OBS_full(k,3:4) = OBS_data(idx,3:4);
    k = k + 1;
end

OBS_full = sortrows(OBS_full,2);

[OBS_dv_n, OBS_DV] = degreeVariance(OBS_full);
OBS_dv_mat = [OBS_dv_n OBS_DV];
save(fullfile(HOME,'Data','DV_OBS.mat'),'OBS_dv_mat');

%%M1 Path to cached M1 SH coefficients
D = 38000;

load(fullfile(HOME,'Data','Model1optdensity','Model1_bounds_iter_05.mat'),'Model1');
M1_V_model_max = segment_2layer_model(Model1.l1.bound, Model1.l2.bound, Model1.l3.bound, Model1.l1.dens, Model1.l2.dens, D, Model1);
[M1_dv_n, M1_DV] = degreeVariance(M1_V_model_max);
M1_dv_matrix = [M1_dv_n M1_DV];
save(fullfile(HOME,'Data','DV_M1'),'M1_dv_matrix')

%% Plot comparison
fig = figure('Color','w');   % Figure background white

semilogy(M1_dv_n(2:end), M1_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance: M1 vs Observed');
legend('M1 Bouguer inversion','Observed','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 12, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_M1_comparison.png'));

%%%M2

load(fullfile(HOME,'Data','M2_V_model.mat'),'M2_V_model');
[M2_dv_n, M2_DV] = degreeVariance(M2_V_model);
M2_dv_matrix = [M2_dv_n M2_DV];
save(fullfile(HOME,'Data','DV_M2'),'M2_dv_matrix')

%% Plot comparison
fig = figure('Color','w');   % Figure background white

semilogy(M2_dv_n(2:end), M2_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance: M2 vs Observed');
legend('M2 Airy','Observed','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 12, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_M2_comparison.png'));

%%%M3

load(fullfile(HOME,'Data','M3_V_Model_Te29000.mat'),'M3_V_Model','T_e');
[M3_dv_n, M3_DV] = degreeVariance(M3_V_Model);
M3_dv_matrix = [M3_dv_n M3_DV];
save(fullfile(HOME,'Data','DV_M3'),'M3_dv_matrix')

load(fullfile(HOME,'Data','Model3opt','M3_V_Model_Te51500.mat'),'M3_V_Model','T_e')
[M3_dv_n_opt, M3_DV_opt] = degreeVariance(M3_V_Model);
M3_dv_matrix_opt = [M3_dv_n_opt M3_DV_opt];
save(fullfile(HOME,'Data','DV_M3_opt'),'M3_dv_matrix_opt')

%% Plot comparison
fig = figure('Color','w');   % Figure background white

semilogy(M3_dv_n(2:end), M3_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance: M3 vs Observed');
legend('M3 Flexure','Observed','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 12, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_M3_comparison.png'));

%% Plot comparison flexure opt
%% Plot comparison flexure opt
fig = figure('Color','w');   % Figure background white

semilogy(M3_dv_n(2:end), M3_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(M3_dv_n_opt(2:end), M3_DV_opt(2:end), 'LineWidth',1.2);
semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance: M3 original and optimised vs Observed');
legend('M3 Flexure','M3 Flexure (optimal)','Observed','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 12, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_M3_opt_comparison.png'));



%% Plot comparison
fig = figure('Color','w');   % Figure background white

semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(M1_dv_n(2:end), M1_DV(2:end), 'LineWidth',1.2);
semilogy(M2_dv_n(2:end), M2_DV(2:end), 'LineWidth',1.2);
semilogy(M3_dv_n(2:end), M3_DV(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance');
legend('Observed', 'M1' , 'M2', 'M3','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 15, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_comparison.png'));

%% Plot comparison
fig = figure('Color','w');   % Figure background white

semilogy(OBS_dv_n(2:end), OBS_DV(2:end), 'LineWidth',1.2);
hold on
semilogy(M1_dv_n(2:end), M1_DV(2:end), 'LineWidth',1.2);
semilogy(M2_dv_n(2:end), M2_DV(2:end), 'LineWidth',1.2);
semilogy(M3_dv_n_opt(2:end), M3_DV_opt(2:end), 'LineWidth',1.2);
hold off

xlim([2 nmax]);
xlabel('Spherical harmonic degree n');
ylabel('Degree Variance (log scale)');
title('Degree Variance');
legend('Observed', 'M1' , 'M2', 'M3 (optimal)','Location','best');

% Force white axes background
ax = gca;
ax.Color = 'w';  % This sets the plotting area background to white

% Force black text and axis lines
set(findall(fig,'Type','text'),'Color','k');
set(findall(fig,'Type','axes'),'XColor','k','YColor','k');
set(findall(fig,'Type','colorbar'),'Color','k');

% Ensure white figure background
set(fig, 'Color', 'w');

% Adjust font size and line width for clarity
set(gca, 'FontSize', 12, 'LineWidth', 1);

saveas(gcf, fullfile(HOME,'Results','DV_comparison_opt.png'));

