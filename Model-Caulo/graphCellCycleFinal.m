%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for plotting results of Caulobacter cell cycle model. 
% Input tout and yout from numerical integration
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

function [] = graphCellCycleFinal(tout, yout)
set(groot, 'DefaultAxesFontSize', 23)
set(groot, 'DefaultFigurePosition', [200 50 600 600])
set(groot, 'DefaultFigureColormap',jet)
set(groot, 'DefaultFigureGraphicsSmoothing','off')
set(groot, 'DefaultAxesLineWidth', 1.2)

% load time and outputs
load("./LoadMatrices/model_result.mat","tout","yout")

%load all experimental data
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
Sup = 18;
DivKp = 19;
I = 20;
II = 21;
III = 22;
tot = 34;

CPLX1 = 23;
CpdR = 24;
CpdRP = 25;
CPLX2 = 26;
RcdA = 27;
CPLX3 = 28;

CtrAP = 29;
CckAP=30; 
cdG=31; 
PleD=32; 
PdeA=33;
PleDP=34;
%data
t =  [5, 30, 60, 90, 120, 150];
dDnaA = [1262	1164 500	496 	1066	1028]./445;
dCcrM = [17	17	24	445	435	101]./445;
dSciP = [2858	451	199	796	2956	2649]./445;
dGcrA = [550	3001	1275	476	725	1493]./445;
dCtrA = [195	612	2251	3216	2486	1617]./445;
% publication 1:
tp1 = [0 20 40 60 80 100 120 140];
pDnaA1 = [4784.305 10537.962 9511.012 3134.012 1447.184 1014.355 3693.083 6339.548];
pGcrA1 = [0 8632.062 16136.355 14185.012 14835.891 8353.184 9923.598 10884.719];
pCtrA1 = [6648.012 0 0 0 4042.598 8955.376 7497.598 6225.962];

%publication 2:
tp2 = [0 20 40 60 80 100 120 140];
pSciP2 = [9239.255 9518.134 4585.305 610.042 161.657 388.435 4424.062 10748.426];
pGcrA2 = [0 444.87 4548.648 6148.669 5947.891 4592.062 5067.406 5809.548];
pCtrA2 = [7177.962 9270.255 1045.406 0 315.506 5709.77 8172.598 7724.355];

%publication 3:
tp3 = [0 16.5 33 49.5 66 82.5 99 115.5 132 148.5];
pCcrM3 = [1314.983 1261.619 2128.518 1517 2385.477 2660.548 2990.276 7318.69 7300.154 4928.518];
pCtrA3 = [4366.669 544.698 678.284 1565.962 3814.74 6857.891 7663.376 13218.74 12383.598 10310.083];

%publication 4:
tp4 = [0 10 20 30 40 50 60 70 80 90 100];
pSciP4 = [8176.062 3082.648 894.456 248.12 248.12 248.12 248.12 855.77 1727.527 2709.941 2193.042];
pCtrA4 = [3465.406 1278.77 928.527 590.991 1960.355 3866.477 4892.77 5509.062 4213.941 4028.77 2851.042];

%publication 5:
tp5 = [0 20 40 60 80 100 120 140];
pDnaA5 = [6925.376 14430.033 13141.619 4829.376 2257.719 2059.497 6289.912 9555.912];
pGcrA5 = [0 8417.062 16478.891 14706.962 16155.255 9149.255 10833.548 12333.669];
pCtrA5 = [8860.134 0 0 0 4023.77 12430.841 10419.669 7582.255];
pCcrM5 = [0 0 0 0 166 3633.548 10993.962 13607.912];

pCpdRP = [9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
tpCpdRP = [20,40,60,80,100,120,140];

pCpdR = [17942.489,19085.196,12028.882,9405.711,9145.711,9438.782,12948.823];
tpCpdR = [0, 20,40,60,80,100,120];

rcda =[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
tpRcdA =[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341];
% CpdR1 Antonio A. Iniesta and Lucy Shapiro. A bacterial control circuit
% integrates polar localization and proteolysis
pCpdRtot = [0.6, 0.41, 0.22, 0.19, 0.19, 0.32, 0.67, 1, 0.83];
tpCpdRtot = [0*0.937, 20*0.937,40*0.937,60*0.937,80*0.937,100*0.937,120*0.937, 140*0.937, 160*0.937];

%PleD and PdeA
tpPdeA = [0 18.75 37.5 56.25 75 93.75 112.5 131.25 150];
dPdeA = [1  0.41  0.158  0.065  0.072  0.173  0.44   0.82  0.998]; %Abel 2011

tpcdG = [0    20    40    60    80   100   120 140];
pcdG = [ 1/3   1   0.55  0.35   0.3  0.25  0.05  0.35]; %Abel 2013

tpcdG2 = [0    20    40    60    80   100   120];
pcdG2 = [0.2292 1.0000 0.5000 0.2833 0.2250 0.1625 0];

% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350));
scaled_dPdeA = (dPdeA - min(dPdeA))/(max(dPdeA)-min(dPdeA))*(max(yout(a:b, PdeA))-min(yout(a:b, PdeA)))+min(yout(a:b, PdeA));
scaled_dCcrM = (dCcrM - min(dCcrM))/(max(dCcrM)-min(dCcrM))*(max(yout(a:b, mCcrM))-min(yout(a:b, mCcrM)))+min(yout(a:b, mCcrM));

% set to 1 to save figures
save_figs = 0;

%% ccrM
figure(1)
line(tout(a:b)-1200, yout(a:b, mCcrM), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-'); %simulated levels
hold on;
box on;
plot(t, scaled_dCcrM, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
xlim([0 150])
ylim([0 1])
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
title('\it{ccrM}', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
yticks([0 0.2 0.4 0.6 0.8 1.0])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/ccrM.eps','Resolution',300)
end

%% dnaA
figure(2)
scaled_dDnaA = (dDnaA - min(dDnaA))/(max(dDnaA)-min(dDnaA))*(max(yout(a:b, mDnaA))-min(yout(a:b, mDnaA)))+min(yout(a:b, mDnaA));
line(tout(a:b)-1200, yout(a:b, mDnaA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-'); % simulated dnaA
hold on;
box on;
plot(t, scaled_dDnaA, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
title('\it{dnaA}', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
xlim([0 150])
ylim([0.5 2.5])
yticks([0.5 1 1.5 2 2.5])
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/dnaA.eps','Resolution',300)
end

%% gcrA
figure(3)
scaled_dGcrA = (dGcrA - min(dGcrA))/(max(dGcrA)-min(dGcrA))*(max(yout(a:b, mGcrA))-min(yout(a:b, mGcrA)))+min(yout(a:b, mGcrA));
line(tout(a:b)-1200, yout(a:b, mGcrA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dGcrA, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
xlim([0 150])
ylim([0 6])
yticks([0 2 4 6])
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
title('\it{gcrA}', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/gcrA.eps','Resolution',300)
end

%% sciP
figure(4)
scaled_dSciP = (dSciP - min(dSciP))/(max(dSciP)-min(dSciP))*(max(yout(a:b, mSciP))-min(yout(a:b, mSciP)))+min(yout(a:b, mSciP));
line(tout(a:b)-1200, yout(a:b, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dSciP, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([0 6])
yticks([0 2 4 6])
title('\it{sciP}', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/sciP.eps','Resolution',300)
end

%% ctrA
figure(5)
scaled_dCtrA = (dCtrA - min(dCtrA))/(max(dCtrA)-min(dCtrA))*(max(yout(a:b, mCtrA))-min(yout(a:b, mCtrA)))+min(yout(a:b, mCtrA));
line(tout(a:b)-1200, yout(a:b, mCtrA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dCtrA, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([0 8])
yticks([0 2 4 6 8])
title('\it{ctrA}', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/ctrA.eps','Resolution',300)
end

%% CcrM
figure(6)
scaled_pCcrM3 = (pCcrM3 - min(pCcrM3))/(max(pCcrM3)-min(pCcrM3))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));
scaled_pCcrM5 = (pCcrM5 - min(pCcrM5))/(max(pCcrM5)-min(pCcrM5))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));

line(tout(a:b)-1200, yout(a:b, CcrM), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tp3, scaled_pCcrM3, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
plot(tp5, scaled_pCcrM5, 'b^', 'MarkerFaceColor', 'b', "MarkerSize", 15)
legend('Sim.','Exp. 1', 'Exp. 2', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
title('CcrM','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
ylim([0 0.8])
yticks([0 0.2 0.4 0.6 0.8])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CcrM.eps','Resolution',300)
end

%% DnaA
figure(7)
scaled_pDnaA1 = (pDnaA1 - min(pDnaA1))/(max(pDnaA1)-min(pDnaA1))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));
scaled_pDnaA5 = (pDnaA5 - min(pDnaA5))/(max(pDnaA5)-min(pDnaA5))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));
line(tout(a:b)-1200, yout(a:b, DnaA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pDnaA1, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
plot(tp5, scaled_pDnaA5, 'b^', 'MarkerFaceColor', 'b', "MarkerSize", 15)
legend('Sim.','Exp. 1', 'Exp. 2', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([0.8 2.6])
yticks([1 1.5 2 2.5])
title('DnaA','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/DnaA.eps','Resolution',300)
end

%% GcrA
figure(8)
scaled_pGcrA1 = (pGcrA1 - min(pGcrA1))/(max(pGcrA1)-min(pGcrA1))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));
scaled_pGcrA2 = (pGcrA2 - min(pGcrA2))/(max(pGcrA2)-min(pGcrA2))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));
line(tout(a:b)-1200, yout(a:b, GcrA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pGcrA1, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
plot(tp2, scaled_pGcrA2, 'b^', 'MarkerFaceColor', 'b', "MarkerSize", 15)
legend('Sim.','Exp. 1', 'Exp. 2', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([3 7])
yticks([3 4 5 6 7])
title('GcrA','FontSize',40)
xlabel('Time (min)','FontSize', 25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/GcrA.eps','Resolution',300)
end

%% SciP
figure(9)
scaled_pSciP2 = (pSciP2 - min(pSciP2))/(max(pSciP2)-min(pSciP2))*(max(yout(a:b, SciP))-min(yout(a:b, SciP)))+min(yout(a:b, SciP));
line(tout(a:b) - 1200, yout(a:b, SciP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tp2, scaled_pSciP2, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
legend('Sim.','Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([1.8 10])
yticks([2 4 6 8 10])
title('SciP','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/SciP.eps','Resolution',300)
end

%% CtrA
figure(10)
CtrAT=yout(a:b, CtrA)+yout(a:b,CtrAP);

scaled_pCtrA3 = (pCtrA3 - min(pCtrA3))/(max(pCtrA3)-min(pCtrA3))*(max(CtrAT)-min(CtrAT))+min(CtrAT);
scaled_pCtrA5 = (pCtrA5 - min(pCtrA5))/(max(pCtrA5)-min(pCtrA5))*(max(CtrAT)-min(CtrAT))+min(CtrAT);
line(tout(a:b)-1200, yout(a:b, CtrA)+yout(a:b, CtrAP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b)-1200,yout(a:b,CtrAP),'Color','k','LineWidth', 3, 'Linestyle', '--');
hold on;
box on;
plot(tp3, scaled_pCtrA3, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
plot(tp5, scaled_pCtrA5, 'b^', 'MarkerFaceColor', 'b', "MarkerSize", 15)
legend('Sim. CtrA','Sim. CtrA~P', 'Exp. CtrA', 'Exp. CtrA', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal');
xlim([0 150])
ylim([0 11])
yticks([0 2 4 6 8 10])
title('CtrA','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CtrA.eps','Resolution',300)
end

%% Hemimethylated_States
figure(11)
box on
title('Hemimethylated States', 'FontSize', 35)
line(tout(a:b) - 1200, yout(a:b, hcori), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b) - 1200, yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b) - 1200, yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');

legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}','location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal');
xlabel('Time (min)','FontSize',25);
ylabel('Probability','FontSize',25);
ylim([0 1])
yticks([0.0 0.2 0.4 0.6 0.8 1.0])
xlim([0 150])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/Hemimethylated_States.eps','Resolution',300)
end

%% DNA synthesis
figure(12)
box on
title('DNA Synthesis', 'FontSize', 40)
line(tout(a:b) - 1200, yout(a:b, Elong), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b) - 1200, yout(a:b, DNA), 'Color', 'm', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b) - 1200, yout(a:b, Count), 'Color', 'b', 'LineWidth', 3, 'Linestyle', '--');

legend( 'Elongation',  'DNA', 'Chromosome','location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal');
xlabel('Time (min)','FontSize',25)
ylabel('Count','FontSize',25)
ylim([0 3])
xlim([0 150])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/DNA_Synthesis.eps','Resolution',300)
end

%% CpdR
figure(13)
scaled_pCpdR = (pCpdR - min(pCpdR))/(max(pCpdR)-min(pCpdR))*(max(yout(a:b, CpdR))-min(yout(a:b, CpdR)))+min(yout(a:b, CpdR));
CpdRtot = yout(a:b, CpdR)+yout(a:b,CpdRP);
scaled_pCpdRtot = (pCpdRtot - min(pCpdRtot))/(max(pCpdRtot)-min(pCpdRtot))*(max(CpdRtot)-min(CpdRtot))+min(CpdRtot);

line(tout(a:b)-1200, yout(a:b,CpdR)+yout(a:b,CpdRP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');

hold on;
box on;
plot(tpCpdRtot, scaled_pCpdRtot, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
legend('Sim.','Exp.','location', 'southoutside','FontSize',19, 'orientation', 'horizontal')
xlim([0 150])
title('CpdR','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CpdR.eps','Resolution',300)
end

%% RcdA
figure(14)
RCDA = yout(:, RcdA);
RCDA = RCDA(a:b);                    % gathering relevant simulated CpdR data
plot(tout(a:b)-1200,  yout(a:b, RcdA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-')
xlim([0 150])
hold on
rcda = (rcda - min(rcda))/(max(rcda)-min(rcda))*(max(RCDA)-min(RCDA))+min(RCDA);
plot(tpRcdA, rcda, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)  % plotting experimental rcda points
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
legend('Sim.', 'Exp.','location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
title('RcdA','FontSize',40)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/RcdA.eps','Resolution',300)
end

%% cdG
figure(15)
scaled_pcdG = (pcdG - min(pcdG))/(max(pcdG)-min(pcdG))*(max(yout(a:b,cdG))-min(yout(a:b,cdG)))+min(yout(a:b,cdG));
scaled_pcdG2 = (pcdG2 - min(pcdG2))/(max(pcdG2)-min(pcdG))*(max(yout(a:b,cdG))-min(yout(a:b,cdG)))+min(yout(a:b,cdG));
line(tout(a:b)-1200, yout(a:b,cdG), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tpcdG, scaled_pcdG, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
hold on;

xlim([0 150])
title('cdG','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
legend('Sim.', 'Exp.','location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/cdG.eps','Resolution',300)
end

%% PleD
figure(16)
tpPleD=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140];
pPleD =[0.097 0.112 0.089 0.067 0.049 0.036 0.027 0.021 0.03 0.057 0.085 0.107 0.123 0.061 0.089];
scaled_pPleD = (pPleD - min(pPleD))/(max(pPleD)-min(pPleD))*(max(yout(a:b,PleD)+yout(a:b,PleDP))-min(yout(a:b,PleD)+yout(a:b,PleDP)))+min(yout(a:b,PleD)+yout(a:b,PleDP));
line(tout(a:b)-1200, yout(a:b, PleD)+yout(a:b, PleDP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
line(tout(a:b)-1200, yout(a:b, PleDP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '--');

hold on;
box on;
plot(tpPleD, scaled_pPleD, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)

legend('Sim. PleD','Sim. PleD~P','Exp.',  'location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
xlim([0 150])
ylim([0.5 1.4])
yticks([0.6 0.8 1 1.2 1.4])
title('PleD','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/PleD.eps','Resolution',300)
end

%% PdeA
figure(17)
line(tout(a:b)-1200, yout(a:b, PdeA), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
hold on;
box on;
plot(tpPdeA, scaled_dPdeA, 'ro', 'MarkerFaceColor', 'r', "MarkerSize", 15)
legend('Sim.','Exp.',  'location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
xlim([0 150])
title('PdeA','FontSize',40)
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/PdeA.eps','Resolution',300)
end

%% CckAP
figure(18)
box on
line(tout(a:b)-1200, yout(a:b, CckAP), 'Color', 'k', 'LineWidth', 3, 'Linestyle', '-');
legend('Sim.', 'location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
title('CckAP','FontSize',40)
xlim([0 150])
xlabel('Time (min)','FontSize',25)
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CckAP.eps','Resolution',300)
end

%% CPLX12
figure(19)
box on
line(tout(a:b)-1200,yout(a:b,CPLX1),'Color', 'k', 'LineWidth', 3, 'Linestyle', '-')
line(tout(a:b)-1200,yout(a:b,CPLX2),'Color', 'b', 'LineWidth', 3, 'Linestyle', '-')
legend('C1','C2', 'location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
xlim([0 150])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CPLX12.eps','Resolution',300)
end

%% CPLX3
figure(20)
box on
line(tout(a:b)-1200,yout(a:b,CPLX3),'Color', 'r', 'LineWidth', 3, 'Linestyle', '-')
legend('C3', 'location', 'southoutside', 'FontSize',19, 'orientation', 'horizontal')
xlim([0 150])
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/CPLX3.eps','Resolution',300)
end

%% DivJ
figure(21)
dataTime = [0 20 40 60 80 100 120 140 ];
dataDivJ=[43316.706 9826.811 14025.225 15523.246 13078.69 12171.447 13689.861 12089.368];
dataDivJ = dataDivJ/max(dataDivJ);
hold on;
box on;
plot(dataTime, dataDivJ, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
legend('Exp.', 'location', 'southoutside', 'FontSize',19,'orientation', 'horizontal')
title('DivJ', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
xlim([0 150])
ylabel('Normalized Conc.','FontSize',25)
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/DivJ.eps','Resolution',300)
end

%% PleC
figure(22)
t_d=1:0.1:150;
 a1 =       80.09;
       b1 =     0.01299;
       c1 =        1.74;
       a2 =       78.77;
       b2 =     0.01333;
       c2 =       4.854;
PleC =  a1*sin(b1*t_d+c1) + a2*sin(b2*t_d+c2);
dataPleC = [84.6 49.2 28.2 42.9 64.1 74.1 80.5 99.1];
dataPleC = dataPleC/max(dataPleC);
dataTime = [0 20 40 80 100 120 140 160];
hold on;
box on;
line(t_d,PleC,'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
hold on
plot(dataTime,dataPleC, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize',15)
title('PleC', 'FontSize', 40)
xlabel('Time (min)','FontSize',25)
xlim([0 150])
ylabel('Normalized Conc.','FontSize',25)
legend('Fitted Curve','Exp.','location', 'southoutside', 'FontSize',19,'orientation', 'horizontal')
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/PleC.eps','Resolution',300)
end

%% bar chart test
figure(23)
f = gcf;  pbaspect([1.3 1 1])
box on;
X = categorical({'DnaA','GcrA','CtrA','CcrM', 'SciP'});
max_DnaA = max(yout(a:b, DnaA));
max_GcrA = max(yout(a:b, GcrA));
% CtrA~P
%max_CtrA = max(yout(a:b, CtrA));
max_CtrA = max(CtrAT);
max_CcrM = max(yout(a:b, CcrM));
max_SciP = max(yout(a:b, SciP));
Y = [max_DnaA, max_GcrA, max_CtrA, max_CcrM, max_SciP];
bar(X,Y,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[1 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/bar_chart.eps','Resolution',300)
end
%% line chart

figure(24)
f = gcf;  pbaspect([1.3 1 1])
set(f, 'DefaultAxesFontSize', 14)

subplot(5,1,1)
avg = ones(1,30);
space1 = zeros(1,90);
space2 = zeros(1,30);
et = 1:1:150;
e1 = [space1,  avg,  space2];
scatter(et, e1, 'sr', 'LineWidth', 3, 'MarkerFaceColor', 'r' );
hold on;
ybar = zeros(size(tout));
threshold = (max(yout(:, CcrM)) + min(yout(:, CcrM)))/2;
ybar(yout(:, CcrM)> threshold )=2;
scatter(tout, ybar, 'sb', 'LineWidth', 3, 'MarkerFaceColor', 'b' );
axis([0 150 0.5 2.5])
grid on 
box on;
set(gca,'xtick',[0:30:150])
set(gca,'xticklabel',[])
set(gca,'ytick',[1,2])
set(gca,'yticklabel',{'Exp.','Sim.'})

yyaxis right
set(gca,'yticklabel',{'', 'CcrM', ''}, 'YColor','black')


subplot(5,1,2)
avg = ones(1,50);
space1 = zeros(1,10);
space2 = zeros(1,90);
e1 = [space1,  avg,  space2];
scatter(et, e1, 'sr', 'LineWidth', 3, 'MarkerFaceColor', 'r' );
hold on;
ybar = zeros(size(tout));
threshold = (max(yout(:, DnaA)) + min(yout(:, DnaA)))/2;
ybar(yout(:, DnaA) > threshold)=2;
scatter(tout, ybar, 'sb', 'LineWidth', 3, 'MarkerFaceColor', 'b' );
axis([0 150 0.5 2.5])
grid on;
box on;
set(gca,'xtick',[0:30:150])
set(gca,'xticklabel',[])
set(gca,'ytick',[1,2])
set(gca,'yticklabel',{'Exp.','Sim.'})

yyaxis right
set(gca,'yticklabel',{'', 'DnaA', ''}, 'YColor','black')

subplot(5,1,3)
avg = ones(1,55);
space1 = zeros(1,45);
space2 = zeros(1,50);
e1 = [space1,  avg,  space2];
scatter(et, e1, 'sr', 'LineWidth', 3, 'MarkerFaceColor', 'r' );
hold on;
ybar = zeros(size(tout));
threshold = (max(yout(:, GcrA)) + min(yout(:, GcrA)))/2;
ybar(yout(:, GcrA)>threshold)=2;
scatter(tout, ybar, 'sb', 'LineWidth', 3, 'MarkerFaceColor', 'b' );
axis([0 150 0.5 2.5])
grid on 
box on;
set(gca,'xtick',[0:30:150])
set(gca,'xticklabel',[])
set(gca,'ytick',[1,2])
set(gca,'yticklabel',{'Exp.','Sim.'})

yyaxis right
set(gca,'yticklabel',{'', 'GcrA', ''}, 'YColor','black')

subplot(5,1,4)
avg = ones(1,15);
avg2 = ones(1,60);
space1 = zeros(1,75);
e1 = [avg, space1, avg2];
scatter(et, e1, 'sr', 'LineWidth', 3, 'MarkerFaceColor', 'r' );
hold on;
ybar = zeros(size(tout));
threshold = (max(yout(:, SciP)) + min(yout(:, SciP)))/2;
ybar(yout(:, SciP)> threshold)=2;
scatter(tout, ybar, 'sb', 'LineWidth', 3, 'MarkerFaceColor', 'b' );
axis([0 150 0.5 2.5])
grid on 
box on;
set(gca,'xtick',[0:30:150])
set(gca,'xticklabel',[])
set(gca,'ytick',[1,2])
set(gca,'yticklabel',{'Exp.','Sim.'})

yyaxis right
set(gca,'yticklabel',{'', 'SciP', ''}, 'YColor','black')

subplot(5,1,5)
avg = ones(1,15);
avg2 = ones(1,90);
space1 = zeros(1,45);
e1 = [avg, space1, avg2];
scatter(et, e1, 'sr', 'LineWidth', 3, 'MarkerFaceColor', 'r' );
hold on;
box on;
ybar = zeros(size(tout));
threshold = (max(yout(:, CtrA)) + min(yout(:, CtrA)))/2;
ybar(yout(:, CtrA)> threshold)=2;
scatter(tout, ybar, 'sb', 'LineWidth', 3, 'MarkerFaceColor', 'b' );
axis([0 150 0.5 2.5])
grid on 
set(gca,'xtick',[0:30:150])
set(gca,'ytick',[1,2])
set(gca,'yticklabel',{'Exp.','Sim.'})

xlabel('Time (min)')

yyaxis right
set(gca,'yticklabel',{'', 'CtrA', ''}, 'YColor','black')

if save_figs == 1
    exportgraphics(f,'./Figures/line_chart.eps','Resolution',300)
end

