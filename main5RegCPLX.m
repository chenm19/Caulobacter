function Results = main5RegCPLX(celltype)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for numericaly simulation of Caulobacter cell cycle model
% This script needs 'events5RegCPLX.m' and 'odes5RegCPLX.m' for running.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% celltype determines which cell cycle (swarmer or stalked) we are tracking.
clc; clf; close all;  

if nargin < 1
   celltype = 0;
end

if celltype == 'stalked'
   celltype = 1;
elseif celltype == 'swarmer'
   celltype = 0;
else
   celltype = 0;
end

% global variable CELLTYPE is used to send the celltype information to
% 'swstmodelwtsimin.m' which is the equations input.
global CELLTYPE;
CELLTYPE = celltype;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%%variable map
global Ini Elong DNA Count hcori hCcrM hCtrA ...
       mCcrM mDnaA mGcrA mSciP mCtrA ...
       CcrM DnaA GcrA SciP CtrA Zring DivKp ...
       I_DnaA I_GcrA I_ccrM tot;
global CPLX1 CpdR CpdRP CPLX2 RcdA CPLX3; 
global CtrAP CckAP  cdG PleD PdeA;
global PleDP; %13/3/2021 RX
global ksZring;

Ini = 1;
Elong = 2;
DNA = 3;
Count = 4;
hcori = 5;
hCcrM = 6;
hCtrA = 7;
mCcrM = 8;
mDnaA = 9;
mGcrA = 10;
mSciP = 11;
mCtrA = 12;
CcrM = 13;
DnaA = 14;
GcrA = 15;
SciP = 16;
CtrA = 17;
Zring = 18;
DivKp = 19;
I_DnaA = 20;
I_GcrA = 21;
I_ccrM = 22;

CPLX1 = 23;
CpdR = 24;
CpdRP = 25;
CPLX2 = 26;
RcdA = 27;
CPLX3 = 28;
%added 1/26/2021 RX
CtrAP = 29;
CckAP=30; 
cdG=31; 
PleD=32; 
PdeA=33;
PleDP=34;
tot = 34;  %3/13/2021 RX

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%&&&&&&&&%
% Initial value of model variables, for a newborn, wild-type swarmer cell
y0 = zeros(tot, 1);
% load  init value
% load('init_4xdnaA.mat','init');
load('init.mat','init');
y0 = init';
y0(Zring) = 0;
ksZring = 0;
% y0(Ini) = 0.02;  

% Load parameter 
param_Complex();
% Copy_of_param_Complex();
% Copy_2_of_param_Complex();
% end of initial values
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Integration parameters
tspan = [0 1500];
tstart = tspan(1);
tfinal = tspan(2);

options = odeset('Events', @events5RegCPLX, 'RelTol', 1e-5, 'AbsTol', 1e-7);
% end of integration parameters
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Numerical Integration
tout = tstart;
yout = y0';
teout = [];
yeout = [];
ieout = [];
tpre = 0;

% loop for continous simulation of the differential equations.
% i is used to counted the interrupt points through the simulation. It can be
% changed to larger to satisfy longer time simulation once more interrupt 
% points are detected.

while tstart < tfinal
	[T, Y, TE, YE, IE] = ode15s(@odes5RegCPLX, [tstart tfinal], y0, options);
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

%end of integration

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% Save the data
%save('swstmodelwtsim30.txt', 'tout', 'yout', '-ASCII');
%Save the initial data at t  = 270 min

save('output.mat', 'tout', 'yout')
%   init = yout(end,:);
%   save('init.mat', 'init')

% mutantPlotter(tout, yout)
% FinalFigGenerator(tout, yout)
graphCellCycle(tout, yout)
%graphCellCycle2(tout, yout)
%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%