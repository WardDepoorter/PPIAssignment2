% main file for the complete GSH circle for a particular model
clear;
close all;
clc;

HOME = pwd;
addpath([HOME '/Data']);
addpath([HOME '/Tools']);

% Model Construction

new_model = 1;

if new_model == 1

  % Construct new model
  
  Model = struct();
  
  Model.number_of_layers = 2;
  Model.name = 'Moon_Test_1';
  
  % Additional variables
  Model.GM = 4.902800118E12; %check this value; from wikipedia
  Model.Re = 1737400;
  Model.geoid = 'none';
  Model.nmax = 179;     
  Model.correct_depth = 0;
  
  % Top layer
  Model.l1.bound = img2mx('Topography/LDEM_4.IMG',4);  % meters with respect to reference sphere
  Model.l1.dens  = 2550;
  
  % Second layer
  Model.l2.bound = -50000+zeros(size(Model.l1.bound));     % meters with respect to reference sphere
  Model.l2.dens  = 3300;	   % Density in kg/m3
  
  % Bottom bound
  Model.l3.bound = -100000;    % meters with respect to reference sphere
  
  % Save model in .mat file for use of the new software
  
  save([HOME '/Data/' Model.name '.mat'],'Model')

else
  % Load previous saved model

  model_name = 'Model';
  load([HOME '/Data/' Model.name '.mat']);
end

%%%%%%%%%%%%%%%%%%% Computation area %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Part that cant be modified %%%%%%%%%%%%%%%%%%%%%%%

latLim =    [-89.5 89.5 1];  % [deg] min latitude, max latitude, resolution latitude (preferable similar to latitude)
lonLim =    [-180 180 1];% [deg] min longitude, max longitude, resolution longitude (preferable similar to latitude)
height =    10.0; % height of computation above spheroid
SHbounds =  [0 179]; % Truncation settings: lower limit, upper limit SH-coefficients used

%%%%%%%%%%%%%% Part that can be modified %%%%%%%%%%%%%%%%%%%%%%%

%% Global Spherical Harmonic Analysis 

tic;
[V_Model] = segment_2layer_model(Model.l1.bound,Model.l2.bound,Model.l3.bound,Model.l1.dens,Model.l2.dens,25000,Model);
toc

%% Global Spherical Harmonic Synthesis

tic;
[data] = model_SH_synthesis(lonLim,latLim,height,SHbounds,V_Model,Model);
toc

%% Save data

save([HOME '/Results/data_' Model.name '_' num2str(SHbounds(1)) '_' num2str(SHbounds(2)) '.mat'],'data','V_Model','Model')
