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

%Te_list = [10e3, 15e3, 20e3, 25e3, 29e3, 35e3, 40e3, 45e3, 46e3, 47e3, 48e3, 49e3, 50e3, 50.5e3, 51e3, 51.5e3, 52e3, 52.5e3, 53e3, 54e3, 55e3, 60e3];
Te_list = [51.5e3];
E = 100e9;
sigma = 0.25;
R = 1738e3;
g = 1.62;

S = load(fullfile(HOME,'Data','DV_OBS.mat'),'OBS_dv_mat');
OBS_DV = S.OBS_dv_mat(:,2);
difference_list = NaN(size(Te_list));

for it = 1:numel(Te_list)

    T_e = Te_list(it);   % current elastic thickness
    fprintf('Running Model 3 for Te = %.0f km\n', T_e/1e3);

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
    
    % File name for this Te
    modelFile = fullfile(HOME, 'Data', 'Model3opt', ...
                         ['M3_V_Model_Te' num2str(T_e) '.mat']);
    
    if isfile(modelFile)
        fprintf('Loading cached M3 model for Te = %.0f km\n', T_e/1e3);
        S = load(modelFile, 'M3_V_Model');
        M3_V_Model = S.M3_V_Model;
    else
        fprintf('Computing M3 model for Te = %.0f km\n', T_e/1e3);
        M3_V_Model = segment_2layer_model( ...
            Model3.l1.bound, Model3.l2.bound, Model3.l3.bound, ...
            Model3.l1.dens, Model3.l2.dens, D, Model3);
    
        save(modelFile, 'M3_V_Model', 'T_e');
    end

    %Degree fitting

    [~, DV_Te_Model] = degreeVariance(M3_V_Model);

    log_obs = log10(OBS_DV(2:nmax));
    log_model = log10(DV_Te_Model(2:nmax));
    valid_idx = isfinite(log_obs) & isfinite(log_model);

    if any(valid_idx)
        difference = sum((log_model(valid_idx) - log_obs(valid_idx)).^2);
    else
        difference = Inf;
    end

    
    difference_list(it) = difference;

    fprintf('T_e = %4.0f km, difference = %.4f\n', T_e/1e3,difference)
    
end

[min_dif, idx_best] = min(difference_list);
Te_best = Te_list(idx_best);
fprintf('Best elastic thickness found: %.0f m\n', Te_best);