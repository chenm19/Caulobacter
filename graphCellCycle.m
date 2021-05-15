%graph cellcycle
function [] = graphCellCycle1(tout, yout)
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
t =  [5, 30, 60, 90, 120, 150]+1200;
dDnaA = [1262	1164 500	496 	1066	1028]./445;
dCcrM = [17	17	24	445	435	101]./445;
dSciP = [2858	451	199	796	2956	2649]./445;
dGcrA = [550	3001	1275	476	725	1493]./445;
dCtrA = [195	612	2251	3216	2486	1617]./445;
% publication 1:
tp1 = [0 20 40 60 80 100 120 140] + 1200;
pDnaA1 = [4784.305 10537.962 9511.012 3134.012 1447.184 1014.355 3693.083 6339.548];
pGcrA1 = [0 8632.062 16136.355 14185.012 14835.891 8353.184 9923.598 10884.719];
pCtrA1 = [6648.012 0 0 0 4042.598 8955.376 7497.598 6225.962];

%publication 2:
tp2 = [0 20 40 60 80 100 120 140] + 1200;
pSciP2 = [9239.255 9518.134 4585.305 610.042 161.657 388.435 4424.062 10748.426];
pGcrA2 = [0 444.87 4548.648 6148.669 5947.891 4592.062 5067.406 5809.548];
pCtrA2 = [7177.962 9270.255 1045.406 0 315.506 5709.77 8172.598 7724.355];

%publication 3:
tp3 = [0 16.5 33 49.5 66 82.5 99 115.5 132 148.5] + 1200;
pCcrM3 = [1314.983 1261.619 2128.518 1517 2385.477 2660.548 2990.276 7318.69 7300.154 4928.518];
pCtrA3 = [4366.669 544.698 678.284 1565.962 3814.74 6857.891 7663.376 13218.74 12383.598 10310.083];

%publication 4:
tp4 = [0 10 20 30 40 50 60 70 80 90 100] + 1200;
pSciP4 = [8176.062 3082.648 894.456 248.12 248.12 248.12 248.12 855.77 1727.527 2709.941 2193.042];
pCtrA4 = [3465.406 1278.77 928.527 590.991 1960.355 3866.477 4892.77 5509.062 4213.941 4028.77 2851.042];

%publication 5:
tp5 = [0 20 40 60 80 100 120 140] + 1200;
pDnaA5 = [6925.376 14430.033 13141.619 4829.376 2257.719 2059.497 6289.912 9555.912];
pGcrA5 = [0 8417.062 16478.891 14706.962 16155.255 9149.255 10833.548 12333.669];
pCtrA5 = [8860.134 0 0 0 4023.77 12430.841 10419.669 7582.255];
pCcrM5 = [0 0 0 0 166 3633.548 10993.962 13607.912];

pCpdRP = [9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
tpCpdRP = [20,40,60,80,100,120,140] + 1200;

pCpdR = [17942.489,19085.196,12028.882,9405.711,9145.711,9438.782,12948.823];
tpCpdR = [0, 20,40,60,80,100,120] + 1200;

rcda =[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
tpRcdA =[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341] + 1200;

%PleD and PdeA
tpPdeA = [0 18.75 37.5 56.25 75 93.75 112.5 131.25 150]+1200;
dPdeA = [1  0.41  0.158  0.065  0.072  0.173  0.44   0.82  0.998]; %Abel 2011

tpcdG = [0    20    40    60    80   100   120 140] + 1200;
pcdG = [ 1/3   1   0.55  0.35   0.3  0.25  0.05  0.35]; %Abel 2013

tpcdG2 = [0    20    40    60    80   100   120] + 1200;
pcdG2 = [0.2292 1.0000 0.5000 0.2833 0.2250 0.1625 0];
% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350));
scaled_dPdeA = (dPdeA - min(dPdeA))/(max(dPdeA)-min(dPdeA))*(max(yout(a:b, PdeA))-min(yout(a:b, PdeA)))+min(yout(a:b, PdeA));
scaled_dCcrM = (dCcrM - min(dCcrM))/(max(dCcrM)-min(dCcrM))*(max(yout(a:b, mCcrM))-min(yout(a:b, mCcrM)))+min(yout(a:b, mCcrM));
figure()
subplot(3,2,1)
line(tout, yout(:, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); %simulated levels
hold on;
box on;
plot(t, scaled_dCcrM, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
% legend('Simulation','Empirical', 'location', 'northwest')
title('\it{ccrM}')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,2)
scaled_dDnaA = (dDnaA - min(dDnaA))/(max(dDnaA)-min(dDnaA))*(max(yout(a:b, mDnaA))-min(yout(a:b, mDnaA)))+min(yout(a:b, mDnaA));
line(tout, yout(:, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); % simulated dnaA
hold on;
box on;
plot(t, scaled_dDnaA, 'ro', 'MarkerFaceColor', 'r')
% legend('Simulation','Empirical');
title('\it{dnaA}')
% xlabel('Time (min)')
xlim([1200 1350])
ylabel('Normalized Concentration')

subplot(3,2,3)
scaled_dGcrA = (dGcrA - min(dGcrA))/(max(dGcrA)-min(dGcrA))*(max(yout(a:b, mGcrA))-min(yout(a:b, mGcrA)))+min(yout(a:b, mGcrA));
line(tout, yout(:, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dGcrA, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
% legend('Simulation','Empirical')
title('\it{gcrA}')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,4)
scaled_dSciP = (dSciP - min(dSciP))/(max(dSciP)-min(dSciP))*(max(yout(a:b, mSciP))-min(yout(a:b, mSciP)))+min(yout(a:b, mSciP));
line(tout, yout(:, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dSciP, 'ro', 'MarkerFaceColor', 'r')
% legend('Simulation','Empirical')
xlim([1200 1350])
title('\it{sciP}')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,5)
scaled_dCtrA = (dCtrA - min(dCtrA))/(max(dCtrA)-min(dCtrA))*(max(yout(a:b, mCtrA))-min(yout(a:b, mCtrA)))+min(yout(a:b, mCtrA));
line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dCtrA, 'ro', 'MarkerFaceColor', 'r')
% legend('Simulation','Empirical', 'location', 'northwest')
xlim([1200 1350])
title('\it{ctrA}')
xlabel('Time (min)')
ylabel('Normalized Concentration')

figure()
scaled_pCcrM3 = (pCcrM3 - min(pCcrM3))/(max(pCcrM3)-min(pCcrM3))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));
scaled_pCcrM5 = (pCcrM5 - min(pCcrM5))/(max(pCcrM5)-min(pCcrM5))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));

subplot(3,2,1)
line(tout, yout(:, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp3, scaled_pCcrM3, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pCcrM5, 'b^', 'MarkerFaceColor', 'b')
% legend('Simulation','Empirical 1', 'Empirical 2')
xlim([1200 1350])
title('CcrM')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,2)
scaled_pDnaA1 = (pDnaA1 - min(pDnaA1))/(max(pDnaA1)-min(pDnaA1))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));
scaled_pDnaA5 = (pDnaA5 - min(pDnaA5))/(max(pDnaA5)-min(pDnaA5))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));
line(tout, yout(:, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pDnaA1, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pDnaA5, 'b^', 'MarkerFaceColor', 'b')
% legend('Simulation','Empirical 1', 'Empirical 2')
xlim([1200 1350])
title('DnaA')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,3)
scaled_pGcrA1 = (pGcrA1 - min(pGcrA1))/(max(pGcrA1)-min(pGcrA1))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));
scaled_pGcrA2 = (pGcrA2 - min(pGcrA2))/(max(pGcrA2)-min(pGcrA2))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));
line(tout, yout(:, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pGcrA1, 'ro', 'MarkerFaceColor', 'r')
plot(tp2, scaled_pGcrA2, 'b^', 'MarkerFaceColor', 'b')
% legend('Simulation','Empirical 1', 'Empirical 2', 'location', 'southeast')
xlim([1200 1350])
title('GcrA')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,4)
scaled_pSciP2 = (pSciP2 - min(pSciP2))/(max(pSciP2)-min(pSciP2))*(max(yout(a:b, SciP))-min(yout(a:b, SciP)))+min(yout(a:b, SciP));
line(tout, yout(:, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp2, scaled_pSciP2, 'ro', 'MarkerFaceColor', 'r')
% legend('Simulation', 'Empirical')
xlim([1200 1350])
title('SciP')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,5)
CtrAT=yout(a:b, CtrA)+yout(a:b,CtrAP);
scaled_pCtrA3 = (pCtrA3 - min(pCtrA3))/(max(pCtrA3)-min(pCtrA3))*(max(CtrAT)-min(CtrAT))+min(CtrAT);
scaled_pCtrA5 = (pCtrA5 - min(pCtrA5))/(max(pCtrA5)-min(pCtrA5))*(max(CtrAT)-min(CtrAT))+min(CtrAT);
line(tout, yout(:, CtrA)+yout(:, CtrAP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
line(tout,yout(:,CtrAP),'Color','k','LineWidth', 2, 'Linestyle', '--');
hold on;
box on;
plot(tp3, scaled_pCtrA3, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pCtrA5, 'b^', 'MarkerFaceColor', 'b')

% legend('Simulation','Empirical 1', 'Empirical 2', 'location', 'northwest')
xlim([1200 1350])
title('CtrA')
xlabel('Time (min)')
ylabel('Normalized Concentration')

figure()
subplot(2,2,1)
line(tout, yout(:, hcori), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p2 = line(tout, yout(:, hCcrM), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, hCtrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
h = legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}', 'Location', 'North');
xlabel('Time (min)')
ylabel('Probability')
xlim([1200 1350])
subplot(2,2,2)
p2 = line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, DNA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p4 = line(tout, yout(:, Count), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout, yout(:, Ini), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
h = legend( 'Elongation',  'DNA', 'Chromosome');
xlabel('Time (min)')
ylabel('Count')
xlim([1200 1350])
subplot(2,2,3)
line(tout, yout(:, hcori), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p2 = line(tout, yout(:, hCcrM), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, hCtrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
h = legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}', 'Location', 'North');
xlim([0 600])
xlabel('Time (min)')
ylabel('Probability')
subplot(2,2,4)
p2 = line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout, yout(:, DNA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p4 = line(tout, yout(:, Count), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
p5 = line(tout, yout(:, Ini), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
h = legend( 'Elongation',  'DNA', 'Chromosome');
ylabel('Count')
xlim([0 600])
ylim([0 2])


figure()
subplot(3,2,1)
scaled_pCpdR = (pCpdR - min(pCpdR))/(max(pCpdR)-min(pCpdR))*(max(yout(a:b, CpdR))-min(yout(a:b, CpdR)))+min(yout(a:b, CpdR));
% scaled_pCpdRP = (pCpdRP - min(pCpdRP))/(max(pCpdRP)-min(pCpdRP))*(max(yout(a:b, CpdRP))-min(yout(a:b, CpdRP)))+min(yout(a:b, CpdRP));
% line(tout, yout(:, CpdRP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
line(tout, yout(:, CpdR), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
% plot(tpCpdRP, scaled_pCpdRP, 'ro', 'MarkerFaceColor', 'r')
plot(tpCpdR, scaled_pCpdR, 'ro', 'MarkerFaceColor', 'r')
% legend('CpdRP','CpdR','Empirical')
xlim([1200 1350])
title('CpdR')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,2)
RCDA = yout(:, RcdA);
RCDA = RCDA(a:b);                    % gathering relevant simulated CpdR data
plot(tout,  yout(:, RcdA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
xlim([1200 1350])
hold on
rcda = (rcda - min(rcda))/(max(rcda)-min(rcda))*(max(RCDA)-min(RCDA))+min(RCDA);
plot(tpRcdA, rcda, 'ro', 'MarkerFaceColor', 'r')  % plotting experimental rcda points
% xlabel('Time (min)')
ylabel('Normalized Concentration')
% legend('Simulation RcdA', 'Empirical')
title('RcdA')

subplot(3,2,3)
scaled_pcdG = (pcdG - min(pcdG))/(max(pcdG)-min(pcdG))*(max(yout(a:b,cdG))-min(yout(a:b,cdG)))+min(yout(a:b,cdG));
scaled_pcdG2 = (pcdG2 - min(pcdG2))/(max(pcdG2)-min(pcdG))*(max(yout(a:b,cdG))-min(yout(a:b,cdG)))+min(yout(a:b,cdG));
line(tout, yout(:,cdG), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpcdG, scaled_pcdG, 'ro', 'MarkerFaceColor', 'r')
hold on
plot(tpcdG2, scaled_pcdG2, 'b*', 'MarkerFaceColor', 'r')
xlim([1200 1350])
title('cdG')
% xlabel('Time (min)')
ylabel('Normalized Concentration')
% legend('Simulation cdG', 'Empirical')
% PleDT=yout(:, PleD)+yout(:, PleDP);
subplot(3,2,4)
tpPleD=[0 10 20 30 40 50 60 70 80 90 100 110 120 130 140]+1200;
pPleD =[0.097 0.112 0.089 0.067 0.049 0.036 0.027 0.021 0.03 0.057 0.085 0.107 0.123 0.061 0.089];
scaled_pPleD = (pPleD - min(pPleD))/(max(pPleD)-min(pPleD))*(max(yout(a:b,PleD)+yout(a:b,PleDP))-min(yout(a:b,PleD)+yout(a:b,PleDP)))+min(yout(a:b,PleD)+yout(a:b,PleDP));
line(tout, yout(:, PleD)+yout(:, PleDP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
line(tout, yout(:, PleDP), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');

hold on;
box on;
plot(tpPleD, scaled_pPleD, 'ro', 'MarkerFaceColor', 'r')

% legend('Simulation PleD','Empirical 1',  'location', 'northwest')
xlim([1200 1350])
title('PleD')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,5)
line(tout, yout(:, PdeA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpPdeA, scaled_dPdeA, 'ro', 'MarkerFaceColor', 'r')
% legend('Simulation PdeA','Empirical 1',  'location', 'northwest')
xlim([1200 1350])
title('PdeA')
% xlabel('Time (min)')
ylabel('Normalized Concentration')

subplot(3,2,6)
line(tout, yout(:, CckAP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% legend('Simulation CckAP',  'location', 'northwest')
title('CckAP')
xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized Concentration')


figure()
subplot(1,2,1)
line(tout,yout(:,CPLX1),'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
line(tout,yout(:,CPLX2),'Color', 'b', 'LineWidth', 2, 'Linestyle', '-')
legend('C1','C2',  'location', 'northwest')
xlim([1200 1350])
subplot(1,2,2)
line(tout,yout(:,CPLX3),'Color', 'r', 'LineWidth', 2, 'Linestyle', '-')
legend('C3',  'location', 'northwest')
xlim([1200 1350])

% figure()
% subplot(1,2,1)
% dataTime = [0 20 40 60 80 100 120 140 ];
% dataDivJ=[43316.706 9826.811 14025.225 15523.246 13078.69 12171.447 13689.861 12089.368];
% dataDivJ = dataDivJ/max(dataDivJ);
% scatter(dataTime,dataDivJ,'red')
%  legend('DivJ data',  'location', 'northwest')
% subplot(1,2,2)
% t_d=1:0.1:150;
%  a1 =       80.09;
%        b1 =     0.01299;
%        c1 =        1.74;
%        a2 =       78.77;
%        b2 =     0.01333;
%        c2 =       4.854;
% PleC =  a1*sin(b1*t_d+c1) + a2*sin(b2*t_d+c2);
% dataPleC = [84.6 49.2 28.2 42.9 64.1 74.1 80.5 99.1];
% dataPleC = dataPleC/max(dataPleC);
% dataTime = [0 20 40 80 100 120 140 160];
%  line(t_d,PleC,'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
%  hold on
%  scatter(dataTime,dataPleC,'red')
%  legend('PleC curve','PleC data',  'location', 'northwest')