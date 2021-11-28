function [] = main5RegCPLX(para)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for numericaly simulation of Caulobacter cell cycle model
% This script needs 'events5RegCPLX.m' and 'odes5RegCPLX.m' for running.
% Input para: 47 optimized parameters 
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clc; clf; close all;  

global dDnaA dCcrM dSciP dGcrA dCtrA pCcrM3 pCcrM5 pDnaA1 pGcrA1 pSciP2 pCtrA3 pCtrA5 ...
        t tp1 tp2 tp3 tp5 tpPleD pPleD tpPdeA dPdeA tpRcdA pCpdR tpCpdR tpcdG pcdG rcda
global Ini Elong DNA Count hcori hCcrM hCtrA ...
       mCcrM mDnaA mGcrA mSciP mCtrA ...
       CcrM DnaA GcrA SciP CtrA Zring DivKp ...
       I_DnaA I_GcrA I_ccrM tot;
global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 
global CtrAP CckAP cdG PleD PdeA;
global PleDP; %13/3/2021 RX
global ksZring;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% Model Setup
% Load other parameters and experimental data 
load_para();
load_data();

% Initial value of model variables
load('init.mat','init');
y0 = init;
y0(Zring) = 0;
ksZring = 0;
% y0(Ini) = 0.02;  

% Integration parameters
tspan = [0 1500];
T = [];
tstart = tspan(1);
tfinal = tspan(2);

options = odeset('Events', @events5RegCPLX, 'RelTol', 1e-5, 'AbsTol', 1e-7);
odefun = @(t, y) odes5RegCPLX(t, y, para);

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% Numerical Integration
tout = tstart;
yout = y0;
teout = [];
yeout = [];
ieout = [];
tpre = 0;
counter = 0;
dt = 0;

% loop for continous simulation of the differential equations.
% i is used to counted the interrupt points through the simulation. It can be
% changed to larger to satisfy longer time simulation once more interrupt 
% points are detected.

while tstart < tfinal
    
	[T, Y, TE, YE, IE] = ode15s(odefun, [tstart tfinal], y0, options);
	% Accumulate output.  
	nt = length(T);

	tout = [tout;T(2:nt)];
	yout = [yout;Y(2:nt, :)];
	
	teout = [teout;TE];
  	yeout = [yeout;YE];
  	ieout = [ieout;IE];
 
  	% Refersh the initial conditions for differential equations once the interrupt 
  	% point is detected.
  	y0 = Y(nt, :);
 
 	if isscalar(IE) == 0
    		IE = 0;
 	end  

 	switch (IE)
 
        case 1
       	y0(Elong) = Y(nt, Elong) + Y(nt, Count)*0.05;
       	y0(DNA) = Y(nt, DNA) + Y(nt, Count)*0.05;
       	y0(Count) = 2; 
       	y0(hcori) = 1;  
        y0(Ini) = 0;
        ksZring = 1/90;
        
        dt = tout(end) - tpre
        tpre = tout(end);
        counter = counter+1;
        
    	case 2
      	y0(hCcrM) = 1;
        
    	case 3
      	y0(hCtrA) = 1;
%         y0(DivKp) = 0;

    	case 4
      	y0(Elong) = 0;

    	case 5
%         y0(Elong) = 0;
      	y0(DNA) = 1;
      	y0(Count) = 1;
      	y0(Zring) = 0;
        ksZring = 0;
%         y0(hcori) = 0; 
%         y0(hCcrM) = 0;
%         y0(hCtrA) = 0;
%         y0(Ini) = 0;
        
        case 6
        y0(DivKp) = 1;
        
%         case 6
%          y0(Ini) = 0.0;
        
	end

	tstart = T(nt);
	if tstart >= tfinal
		break;
	end

end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% Save and plot the results

save('model_result.mat', 'tout', 'yout')

%  init = yout(end,:);
%  save('init.mat', 'init')

%  mutantPlotter(tout, yout)
%  FinalFigGenerator(tout, yout)
graphCellCycle(tout, yout)

%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%