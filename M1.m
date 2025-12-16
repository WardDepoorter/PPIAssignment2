%M1: run 2_layermodel first to set initial model setup
clear;
close all;
clc;

load(["Data/Moon_Test_2.mat"]);%make sure to run 2_layer file first with the desired params

% x0 = [2800, 38000];
% options = optimset('Display', 'iter');
% result = fminsearch(@M1_MSE_calculator,x0, options );
% 
x = zeros(1,2);

t_crust = 38000;
x(2) = t_crust;
result = []; % Initialize result as an empty array
i = 1;
for rho = 2500:100:3500
    x(1) = rho;
    result(i) = M1_MSE_calculator(x); % Store the result of each iteration
    i = i + 1;
end
%%plot result array
figure;
plot(2500:100:3500, result, '-o');
xlabel('Density (kg/m^3)');
ylabel('Mean Squared Error');
title('MSE vs Density');
grid on;