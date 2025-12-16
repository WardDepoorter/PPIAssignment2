%plots grail datafile; see
%https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/rsdmap/gggrx_0660pm_anom_l320.lbl
%for info. contains gravity anomaly in miligal 

% filename='Data/gggrx_0660pm_anom_l320.img'
% fid = fopen(filename, 'r', 'ieee-le');
% dataa = fread(fid, [1440, 721], 'single=>single'); 
% fclose(fid);
% Mx = double(dataa');
% disp(['Min: ', num2str(min(Mx(:))), ', Max: ', num2str(max(Mx(:)))]) 
% lowres = compress_to_1deg(Mx, 4)
% 
% 
% load('Results/data_Moon_Test_1_0_179.mat')
% 
% lon = data.grd.lon(1,:);
% lats = data.grd.lat(:,1);
% 
% error = lowres - flipud(((data.vec.R)).*1e5);
% 



% Plot
% figure;
% imagesc([-180 180], -[-90 90], error);
% set(gca, 'YDir', 'normal');  % so latitude increases upward
% demcmap(error, 128); % Using a better colormap for topography visualization
% colorbar();
% 
% %title('Crustal Thickness');
% xlabel('Longitude (deg)');
% ylabel('Latitude (deg)');
%1738,4902.79980693169,7.74304189736151e-06,660,660,1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function synthesizes the SH from Grail to plot and compare.

clear;
close all;
clc;

mx = txt2mx('GRAIL/Grail_sh.txt');
mx = mx(:, 1:end-2); % Remove the last two columns containing uncertainties
V_Model = mx;
Model = struct();
  Model.name = 'Graildata';
  % Additional variables
  Model.GM = 4.9028E12;
  Model.Re = 1738*1000;
  Model.nmax = 330;


latLim =    [-89.5 89.5 0.25];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 0.25];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    0; % height of computation above spheroid
SHbounds =  [0 330]; % Truncation settings: lower limit, upper limit SH-coefficients used

[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);

save(['Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'],'data')
