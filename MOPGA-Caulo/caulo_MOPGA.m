%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for multiobjective optimizations for Caulobacter cell cycle model. 
% 47 parameters need to be optimized
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clc; clear; close all;
warning('off','MATLAB:ode15s:IntegrationTolNotMet');

% Load all parameters and experimental data 
load_para();
load_data();

% load starting point
para = [1.4 1.4 1.4 0.1 0.0667 0.256 0.08 0.242 0.06 5.6 ...
         0.6 0.5 0.04 0.99 0.09 0.083 0.065 0.028 0.085 0.118 ...
         0.06 0.0432 0.06 0.6 3 0.7 1.5 1 1 1.1 ...
         1 0.15 0.2 140 2 0.01 1 0.1 0.15 0.04 ...
         0.04 0.01 0.5 5 1 0.1 1];

ub = para.*4;
lb = para./4;

% load good starting point
% load('best4S.mat')
% para = val(2,:);


% Multiobjective GA
options = optimoptions('gamultiobj','MaxGenerations', 30,'PopulationSize', 50, ...
                       'InitialPopulation', para, 'PlotFcn',@gaplotpareto);

[val, fval] = gamultiobj(@fitness, 47,[],[],[],[],lb,ub,options)

save("./LoadMatrices/opt_MultiGA_Outputs.mat", "val", "fval");

