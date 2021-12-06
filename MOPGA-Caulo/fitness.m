%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for objective functions of parameter optimizations in Caulobacter 
% cell cycle model. 
% This script needs 'events5RegCPLX.m' and 'odes5RegCPLX.m' for running.
% Input para: 47 parameters need to be optimized
% Output error: two objective function values
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

function error = fitness(para)

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

% Load initial value of model variables
load('init.mat','init'); 
y0 = init;
y0(Zring) = 0;
ksZring = 0;

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
        
        dt = tout(end) - tpre; %time between integrations
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
        
        case 6
        y0(DivKp) = 1;
%         
	end

	tstart = T(nt);
	if tstart >= tfinal
		break;
	end

end

%% Calculate the MSE between simulated results and experimental data

% find indices of tspan corresponding to scaled experimental times
% the eighth cell cycle
t_ = interp1(tout,1:length(tout), t+1200,'nearest');  
tp1_ = interp1(tout,1:length(tout), tp1+1200,'nearest');
tp2_ = interp1(tout,1:length(tout), tp2+1200,'nearest');
tp3_ = interp1(tout,1:length(tout), tp3+1200,'nearest'); 
tp5_ = interp1(tout,1:length(tout), tp5+1200,'nearest');
tpPleD_ = interp1(tout,1:length(tout), tpPleD+1200,'nearest');
tpPdeA_ = interp1(tout,1:length(tout), tpPdeA+1200,'nearest');
tpRcdA_ = interp1(tout,1:length(tout), tpRcdA+1200,'nearest');
tpCpdR_ = interp1(tout,1:length(tout), tpCpdR+1200,'nearest');
tpcdG_ = interp1(tout,1:length(tout), tpcdG+1200,'nearest');

%--------------mRNA-------------
scaled_dGcrA = (dGcrA - min(dGcrA))/(max(dGcrA)-min(dGcrA))*(6.7-1)+1;
mGcrA_err = mse(yout(t_, mGcrA) - scaled_dGcrA');
scaled_dDnaA = (dDnaA - min(dDnaA))/(max(dDnaA)-min(dDnaA))*(2.8-1)+1;
dDnaA_err = mse(yout(t_, mDnaA) - scaled_dDnaA');

scaled_dCcrM = (dCcrM - min(dCcrM))/(max(dCcrM)-min(dCcrM));
dCcrM_err = mse(yout(t_, mCcrM) - scaled_dCcrM');

scaled_dCtrA = (dCtrA - min(dCtrA))/(max(dCtrA)-min(dCtrA))*(7.2-.4)+.4;
dCtrA_err = mse(yout(t_, mCtrA) - scaled_dCtrA');

scaled_dSciP = (dSciP - min(dSciP))/(max(dSciP)-min(dSciP))*(6.6-.4)+.4;
dSciP_err = mse(yout(t_, mSciP) - scaled_dSciP');

%--------------Proteins-------------
scaled_pCcrM3 = (pCcrM3 - min(pCcrM3))/(max(pCcrM3)-min(pCcrM3));
CcrM_err = mse(yout(tp3_, CcrM) - scaled_pCcrM3');

% scaled_pCcrM5 = (pCcrM5 - min(pCcrM5))/(max(pCcrM5)-min(pCcrM5));
% CcrM_err = mse(yout(tp5_, CcrM) - scaled_pCcrM5');

scaled_pDnaA1 = (pDnaA1 - min(pDnaA1))/(max(pDnaA1)-min(pDnaA1))*(2);
DnaA_err = mse(yout(tp1_, DnaA) - scaled_pDnaA1');

scaled_pGcrA1 = (pGcrA1 - min(pGcrA1))/(max(pGcrA1)-min(pGcrA1))*(6);
GcrA_err = mse(yout(tp1_, GcrA) - scaled_pGcrA1');

scaled_pSciP2 = (pSciP2 - min(pSciP2))/(max(pSciP2)-min(pSciP2))*(12);
SciP_err = mse(yout(tp2_, SciP) - scaled_pSciP2');

% scaled_pCtrA3 = (pCtrA3 - min(pCtrA3))/(max(pCtrA3)-min(pCtrA3))*(10);
% ctrA = yout(tp3_, CtrA)+ yout(tp3_, CtrAP);
% CtrA_err = mse(ctrA - scaled_pCtrA3');

scaled_pCtrA5 = (pCtrA5 - min(pCtrA5))/(max(pCtrA5)-min(pCtrA5))*(10);
ctrA = yout(tp5_, CtrA)+ yout(tp5_, CtrAP);
CtrA_err = mse(ctrA - scaled_pCtrA5');

%---------------Protease-------------------
PleD_err = mse(yout(tpPleD_, PleD) - pPleD);
PdeA_err = mse(yout(tpPdeA_, PdeA) - dPdeA);
RcdA_err = mse(yout(tpRcdA_, RcdA) - rcda);

scaled_pCpdR = (pCpdR - min(pCpdR))/(max(pCpdR)-min(pCpdR))*(1.2);
CpdR_err = mse(yout(tpCpdR_, CpdR) - scaled_pCpdR);
cdG_err = mse(yout(tpcdG_, cdG) - pcdG);

alpha = 1;
beta = 1;
theta = 1;

mRNA_err = alpha*(mGcrA_err + dCcrM_err + dDnaA_err + dCtrA_err + dSciP_err);
protein_err = beta*(CcrM_err + DnaA_err + GcrA_err + SciP_err + CtrA_err);
protease_err = theta*(PleD_err + PdeA_err + RcdA_err + CpdR_err + cdG_err);

error = zeros(2, 1); 
error(1) = (mRNA_err + protein_err + protease_err)/15; 
error(2) = abs(dt-150);

% if abs(dt-150) < 2
%     error(2) = 0;
% else
%     error(2) = abs(dt-150);
% end

end
