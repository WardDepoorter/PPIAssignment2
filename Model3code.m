%% Setup

clear;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

fprintf('Loading models\n')


%% Load  Moon data

%LOLA topography radius 1737400 m
LOLAtop = 'LDEM_4.IMG';
resolution_lola = 4;
f = fopen(LOLAtop,"r","ieee-le");
r_moon = 0.5 * fread(f,[360*resolution_lola Inf],'int16')' + 1737400;
disp(size(r_moon))
fclose(f);
 
%Grail geoid radius 1378000 m gggrx_1200a_geoid_l660.lbl
GRAILgeoid = 'GRAILgeoid.img';
resolution_GRAIL = 16;
f = fopen(GRAILgeoid,"r","ieee-le");
geoid = fread(f, [360*resolution_lola Inf],'float')' + 1738000;
fclose(f);

%Grail SHA data
%% Gravity data Grail gggrx_0900c_sha.tab
data = readmatrix('GRAILdata2.txt');

[rows, cols] = size(r_moon);

target_res = 1;

%Check if dimensions of the moon are okay with the target resoluation
if mod(rows, target_res) ~= 0 || mod(cols, target_res) ~= 0
    error('Input dimensions are not compatible with the target resolution')
end

%Make the smaller matrix

rows_small = rows / target_res;
cols_small = cols / target_res;
geoid_small = zeros(rows_small, cols_small);

for i = 1:rows_small
    for j = 1:cols_small
        row_index = (i-1)*target_res + 1 : i*target_res;
        col_index = (j-1)*target_res + 1 : j*target_res;
        part = geoid(row_index, col_index);
        geoid_small(i, j) = mean(part(:));
    end
end

fprintf('Starting calculations\n')

topography = r_moon - geoid_small;
nmax = 180;

%remove the header of SHA data
data(1,:) = [];

V = data(:,1:4);
V = [0 0 1 0; V];
Vmax = (nmax + 1)*(nmax + 2)/2;
V = V(1:Vmax, :);
V = sortrows(V,2);

%% Convert to spherical harmonic coefficient matrix
lmax = max(V(:,1));
field = zeros(lmax+1, 2*lmax+1);

for k = 1:size(V,1)
    n = V(k,1);
    m = V(k,2);
    C = V(k,3);
    S = V(k,4);

    field(n+1, lmax+1 + m) = C;
    if m ~= 0
        field(n+1, lmax+1 - m) = S;
    end
end

field(1,:) = 0; %C0,0 is zero

fprintf('Completed calculation of SC\n')

fprintf('Model is being made now\n')

%Crust
D = 38000;
D_topo = -D * ones(size(topography));
D_moon = -D * ones(size(r_moon));

rho_c = 2800;
rho_m = 3300;

r_airy = -topography * rho_c / (rho_m - rho_c);

% Cache file for Airy root
cacheFile = fullfile(HOME, 'Data', 'r_airy.mat');

if isfile(cacheFile)
    fprintf('Loading cached Airy root...\n');
    K = load(cacheFile, 'r_airy');
    r_airy = K.r_airy;
else
    fprintf('Saving Airy root...\n');
    save(cacheFile, 'r_airy', '-v7.3');
    fprintf('Airy root saved.\n');
end

fprintf('Flexural respons calculations')

cs = GSHA(r_airy,nmax);
sc = cs2sc(cs);

n = 1:size(sc,1);
E = 100e9;
T_e = 29e3;
sigma = 0.25;
R = 1738e3;
g = 1.62;

% Flexural rigidity
D_flex = E * T_e^3 / (12 * (1 - sigma^2));

% Degree vector (must be column)
n = (0:size(sc,1)-1)';   % e.g. 0:180 → 181×1

% Wavenumber (element-wise!)
k = (2*n + 1) ./ (2*R);

% Flexural response function (degree-dependent)
phi = (1 + (D_flex ./ ((rho_m - rho_c) * g)) .* k.^4).^(-1);

% Apply flexural filter to spherical harmonics
sc_phi = zeros(size(sc));

for i = 1:size(sc,2)
    sc_phi(:,i) = sc(:,i) .* phi;
end

resolution = 16;
bound_dist = 1/(2*resolution);

latlim = [(-90 + bound_dist), (90 - bound_dist), (1/resolution)];
lonlim = [(0  + bound_dist), (360 - bound_dist), (1/resolution)];

[nlat, nlon] = size(topography);

lonT = linspace(0 + 1/(2*resolution), 360 - 1/(2*resolution), nlon);
latT = linspace(-90 + 1/(2*resolution), 90 - 1/(2*resolution), nlat);
height = 0;
SHbounds = [2 nmax];

gshs_flex = GSHS(sc_phi, lonT, 90-latT, nmax);

Model3 = struct();
Model3.number_of_layers = 2;
Model3.name = 'Two_layered_moon_crust_model3';
Model3.GM = 4.90280011526323e12;
Model3.Re = 1738000;
Model3.geoid = 'none';
Model3.nmax = nmax;
Model3.l1.bound = topography;
Model3.l1.dens = rho_c;
Model3.l2.bound = gshs_flex - D*ones(size(topography));
Model3.l2.dens = rho_m;
Model3.l3.bound = -200000;

fprintf('Model has been made\n')


% Cache file for observed gravity synthesis
cacheFile = fullfile(HOME, 'Data', 'model_data_obs_3.mat');

if isfile(cacheFile)
    fprintf('Loading cached observed gravity model...\n');
    S = load(cacheFile, 'model_data_3');
    model_data_3 = S.model_data_3;
else
    fprintf('Computing observed gravity model...\n');
    model_data_3 = model_SH_synthesis(lonlim, latlim, height, SHbounds, V, Model3);
    save(cacheFile, 'model_data_3', '-v7.3');
    fprintf('Observed gravity model saved.\n');
end

gravity_data_3 = flipud(model_data_3.vec.R) * 1e5; %unit is mGal

fprintf('Starting modeled GSH calculations\n')

% Build cache filename with key parameters
cacheFile_model = fullfile(HOME, 'Data', sprintf('V_model_data_3_nmax%d_D%d_l3_%d_rhoc%d_rhom%d.mat', nmax, D, abs(Model3.l3.bound), rho_c, rho_m));
if isfile(cacheFile_model)
    fprintf('Loading cached modeled gravity data...\n');
    S = load(cacheFile_model, 'V_model_data_3');
    V_model_data_3 = S.V_model_data_3;
else
    fprintf('Computing modeled gravity data...\n');

    % Build SH coefficients for the layered model
    V_model_3 = segment_2layer_model(Model3.l1.bound, Model3.l2.bound, Model3.l3.bound, Model3.l1.dens, Model3.l2.dens, D, Model3);%%

    % Synthesize gravity field
    V_model_data_3 = model_SH_synthesis(lonlim, latlim, height, SHbounds, V_model_3, Model3);%

    % Save to cache%
    save(cacheFile_model, 'V_model_data_3', '-v7.3');
    fprintf('Modeled gravity data saved.\n');
end

gravity_data_model_3 = flipud(V_model_data_3.vec.R) * 1e5; % mGal
lon = V_model_data_3.grd.lon(1,:);
lat = V_model_data_3.grd.lat(:,1);

fprintf('calculating the difference\n')
gravity_dif_3 = gravity_data_3 - gravity_data_model_3;
disp([max(abs(gravity_dif_3(:))), min(gravity_dif_3(:))]);


fprintf('plotting now\n')

figure;
ax1 = subplot(2,2,1);
imagesc(lon,-lat,gravity_data_model_3);
c=colorbar;
colormap(ax1, 'parula')
hold on
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Modeled gravity'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

ax2 = subplot(2,2,2);
imagesc(lon,-lat,gravity_data_3);
c=colorbar;
colormap(ax2, 'parula')
hold on
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Observed gravity'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

ax3 = subplot(2,2,3);
imagesc(lon,-lat,gravity_dif_3);
c=colorbar;
colormap(ax3, 'parula')
hold on
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Residual gravity'])
ylabel(c,'mGal')
set(gca,'YDir','normal')

thickness = topography + D*ones(size(topography)) - gshs_flex;

ax4 = subplot(2,2,4);
imagesc(lon,-lat,thickness/1e3);
c=colorbar;
colormap(ax4, 'turbo')
hold on
xlim([min(lon) max(lon)])
ylim([min(lat) max(lat)])
hold off
xlabel('Longitude [^o]')
ylabel('Latitude [^o]')
title(['Crust thickness'])
ylabel(c,'km')
set(gca,'YDir','normal')

t = sgtitle('Model 3');
t.Color = 'k';

% Force all text and axes to black (for readable saved figures)
set(findall(gcf, 'Type', 'text'), 'Color', 'k');
set(findall(gcf, 'Type', 'axes'), 'XColor', 'k', 'YColor', 'k');
set(findall(gcf, 'Type', 'colorbar'), 'Color', 'k');

% White background
set(gcf, 'Color', 'w');

% Cache file for Model1 figure
figCacheFile = fullfile(HOME, 'Results', 'Model3.png');

if isfile(figCacheFile)
    fprintf('Model3 figure already exists. Skipping save.\n');
else
    fprintf('Saving Model3 figure...\n');
    set(gcf, 'Color', 'w');      % white background
    print(gcf, figCacheFile, '-dpng', '-r300');
    fprintf('Model3 figure saved.\n');
end

