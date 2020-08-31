%% Ploting the time courses of model variables

set(groot, 'DefaultAxesFontSize', 20)
%set(groot,'Defaultfigureposition',[500 500 500 500])

set(groot,'DefaultFigureColormap',jet)
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                                   .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureGraphicsSmoothing','off')
box on

close all; clf;

load('output.mat')

save_figs = 0;  %switch for saving: 1 = save; 0 = no save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Five transcriptional regulators
% figure(11)
% title('Four transcriptional regulators from Shenghua model')
% line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% p2 = line(tout, yout(:, DnaA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
% p3 = line(tout, yout(:, GcrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
% p4 = line(tout, yout(:, CcrM), 'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
% p5 = line(tout, yout(:, SciP), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
% h = legend('CtrA', 'DnaA',  'GcrA', 'CcrM', 'SciP', 'Location', 'North');

% figure(22)
% title('Five transcriptional regulators from Shenghua model')
% p1 = line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% p2 = line(tout, yout(:, mDnaA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
% p3 = line(tout, yout(:, mGcrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
% p4 = line(tout, yout(:, mCcrM), 'Color', 'y', 'LineWidth', 2, 'Linestyle', '-');
% p5 = line(tout, yout(:, mSciP), 'Color', 'g', 'LineWidth', 2, 'Linestyle', '-');
% h = legend('mCtrA', 'mDnaA',  'mGcrA', 'mCcrM', 'mSciP','Location', 'North');

% hold on;

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
tot = 28;

CPLX1 = 23;
CpdR = 24;
CpdRP = 25;
CPLX2 = 26;
RcdA = 27;
CPLX3 = 28;

t =  [5, 30, 60, 90, 120, 150] + 1200;
dDnaA = [1262	1164 500	496 	1066	1028]./445;
dCcrM = [17	17	24	445	435	101]./445;
dSciP = [2858	451	199	796	2956	2649]./445;
dGcrA = [550	3001	1275	476	725	1493]./445;
dCtrA = [195	612	2251	3216	2486	1617]./445;

%% old protein data
tp_old = [0, 20, 40, 60, 80, 100, 120, 140] + 1200;
tpCcrM_old = [0, 20, 40, 60, 80, 100, 120] + 1200;
tp2_old = [9, 27, 44, 62, 80, 98, 117, 134, 152] + 1200;

pDnaA_old = 2*[0.4567	0.9857	0.9147	0.3154	0.128	0.0793	0.4203	0.6582];
pGcrA_old = 6*[0	0.5021	1	0.884	1	0.554	0.662	0.7388];
pCtrA1_old = 7*[0.7791	0	0	0	0.4682	1	0.8591	0.6672];
pCcrM_old = [0.979012341	0.245634747	0.105451343	0.054014407	0.031971569	0.078998499	1];
pSciP_old = 0.7*[7.757911401	7.931725156	4.139518022	1	0.446070124	0.429447455	0.761125117	3.856508369	8.881823903];
pCtrA2_old = [0.8158	0.614	0.1881	0.0323	0.3822	0.7033	0.8924	1	0.782];
pCpdRP = [9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
tpCpdRP = [20,40,60,80,100,120,140] + 1200;

rcda =[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
tpRcdA =[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341] + 1200;

%% new protein data:

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

% Divkp data
T = 150;        % period of Caulobacter
t_d = rem(tout,T); % return remainder after division t/T
DivKP = 0.004218 * t_d + 0.5779;
pDivKP = [0.6311 0.7972 0.8103 0.916 1];
tpDivKP = [20 40 60 80 100] + 1200;

% cdG data:
cdG= 0.3233 * sin(pi/75 * t_d + 0.3469) + 0.3363;
cdG(cdG<0) = 0;
tpcdG = [0    20    40    60    80   100   120] + 1200;
pcdG = [0.2292 1.0000 0.5000 0.2833 0.2250 0.1625 0];

%%
% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350));

%% Plotting concentration of 5 main regulators mRNA against time

% ccrM mRNA plot
% scaling data:
scaled_dCcrM = (dCcrM - min(dCcrM))/(max(dCcrM)-min(dCcrM))*(max(yout(a:b, mCcrM))-min(yout(a:b, mCcrM)))+min(yout(a:b, mCcrM));

figure(1)
 
line(tout(a:b), yout(a:b, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); %simulated levels
hold on;
box on;
plot(t, scaled_dCcrM, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
legend('Simulation','Empiracle', 'location', 'northwest')
title('\it{ccrM}')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ccrM.eps','Resolution',300)
end

% ccrM mRNA full plot
% figure()
%  
% line(tout, yout(:, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(t, dCcrM, 'ro', 'MarkerFaceColor', 'r')
% plot(t, scaled_dCcrM, 'b^', 'MarkerFaceColor', 'b')
% title('\it{ccrM}')
% xlabel('Time (min)')
% ylabel('Scaled ccrM')
% f = gcf;
% xlim([0 1500])
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/ccrM_full.eps','Resolution',300)
% end

% dnaA mRNA plot
% scaled data:
scaled_dDnaA = (dDnaA - min(dDnaA))/(max(dDnaA)-min(dDnaA))*(max(yout(a:b, mDnaA))-min(yout(a:b, mDnaA)))+min(yout(a:b, mDnaA));

figure()
 
line(tout(a:b), yout(a:b, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); % simulated dnaA
hold on;
box on;
plot(t, scaled_dDnaA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Empiracle');
title('\it{dnaA}')
xlabel('Time (min)')
xlim([1200 1350])
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/dnaA.eps','Resolution',300)
end

% plot dnaA full time plot
% figure()
%  
% line(tout, yout(:, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(t, dDnaA, 'ro', 'MarkerFaceColor', 'r')
% plot(t, scaled_dDnaA, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('\it{dnaA}')
% xlabel('Time (min)')
% xlim([0 1500])
% ylabel('Scaled dnaA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/dnaA_full.eps','Resolution',300)
% end


% gcra mRNA plot
% scaled data:
scaled_dGcrA = (dGcrA - min(dGcrA))/(max(dGcrA)-min(dGcrA))*(max(yout(a:b, mGcrA))-min(yout(a:b, mGcrA)))+min(yout(a:b, mGcrA));

figure()
 
line(tout(a:b), yout(a:b, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dGcrA, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
legend('Simulation','Empiracle')
title('\it{gcrA}')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/gcrA.eps','Resolution',300)
end

% gcra mRNA full time subplot
% figure()
%  
% line(tout, yout(:, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(t, dGcrA, 'ro', 'MarkerFaceColor', 'r')
% plot(t, scaled_dGcrA, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('\it{gcrA}')
% xlabel('Time (min)')
% ylabel('Scaled gcrA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/gcrA_full.eps','Resolution',300)
% end

% sciP mRNA plot
% scaled data:
scaled_dSciP = (dSciP - min(dSciP))/(max(dSciP)-min(dSciP))*(max(yout(a:b, mSciP))-min(yout(a:b, mSciP)))+min(yout(a:b, mSciP));

figure()

line(tout(a:b), yout(a:b, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dSciP, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Empiracle')
xlim([1200 1350])
title('\it{sciP}')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/sciP.eps','Resolution',300)
end

% sciP mRNA full time subplot
% figure()
%  
% line(tout, yout(:, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(t, dSciP, 'ro', 'MarkerFaceColor', 'r')
% plot(t, scaled_dSciP, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('\it{sciP}')
% xlabel('Time (min)')
% ylabel('Scaled sciP')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/sciP_full.eps','Resolution',300)
% end

% ctrA mRNA plot
%scaled data:
scaled_dCtrA = (dCtrA - min(dCtrA))/(max(dCtrA)-min(dCtrA))*(max(yout(a:b, mCtrA))-min(yout(a:b, mCtrA)))+min(yout(a:b, mCtrA));

figure()
 
line(tout(a:b), yout(a:b, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, scaled_dCtrA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Empiracle', 'location', 'northwest')
xlim([1200 1350])
title('\it{ctrA}')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ctrA.eps','Resolution',300)
end

% ctra mRNA full time plot
% figure()
%  
% line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(t, dCtrA, 'ro', 'MarkerFaceColor', 'r')
% plot(t, scaled_dCtrA, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('\it{ctrA}')
% xlabel('Time (min)')
% ylabel('Scaled crrA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/ctrA_full.eps','Resolution',300)
% end

%% Plotting concentration of 5 main regulator proteins (and protease complex) against time
% Ccrm protein plot
scaled_pCcrM3 = (pCcrM3 - min(pCcrM3))/(max(pCcrM3)-min(pCcrM3))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));
scaled_pCcrM5 = (pCcrM5 - min(pCcrM5))/(max(pCcrM5)-min(pCcrM5))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM));

figure()
 
line(tout(a:b), yout(a:b, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp3, scaled_pCcrM3, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pCcrM5, 'b^', 'MarkerFaceColor', 'b')
legend('Simulation','Empiracle Source 1', 'Empiracle Source 2')
xlim([1200 1350])
title('CcrM')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Ccrm_protein.eps','Resolution',300)
end

% Ccrm protein full time plot
% figure()
%  
% line(tout, yout(:, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tp3, scaled_pCcrM3, 'ro', 'MarkerFaceColor', 'r')
% plot(tp5, scaled_pCcrM5, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('CcrM')
% xlabel('Time (min)')
% ylabel('Scaled CcrM')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/Ccrm_protein_full.eps','Resolution',300)
% end

% DnaA protein plot
scaled_pDnaA1 = (pDnaA1 - min(pDnaA1))/(max(pDnaA1)-min(pDnaA1))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));
scaled_pDnaA5 = (pDnaA5 - min(pDnaA5))/(max(pDnaA5)-min(pDnaA5))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA));

figure()
 
line(tout(a:b), yout(a:b, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pDnaA1, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pDnaA5, 'b^', 'MarkerFaceColor', 'b')
legend('Simulation','Empiracle Source 1', 'Empiracle Source 2')
xlim([1200 1350])
title('DnaA')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DnaA_protein.eps','Resolution',300)
end

% DnaA protein full time plot
% figure()
%  
% line(tout, yout(:, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tp1, scaled_pDnaA1, 'ro', 'MarkerFaceColor', 'r')
% plot(tp5, scaled_pDnaA5, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('DnaA')
% xlabel('Time (min)')
% ylabel('Scaled DnaA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/DnaA_protein_full.eps','Resolution',300)
% end

% GcrA protein plot
scaled_pGcrA1 = (pGcrA1 - min(pGcrA1))/(max(pGcrA1)-min(pGcrA1))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));
scaled_pGcrA2 = (pGcrA2 - min(pGcrA2))/(max(pGcrA2)-min(pGcrA2))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA));

figure()
 
line(tout(a:b), yout(a:b, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp1, scaled_pGcrA1, 'ro', 'MarkerFaceColor', 'r')
plot(tp2, scaled_pGcrA2, 'b^', 'MarkerFaceColor', 'b')
legend('Simulation','Empiracle Source 1', 'Empiracle Source 2', 'location', 'southeast')
xlim([1200 1350])
title('GcrA')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/GcrA_protein.eps','Resolution',300)
end

%GcrA protein full time subplot
% figure()
%  
% line(tout, yout(:, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tp1, scaled_pGcrA1, 'ro', 'MarkerFaceColor', 'r')
% plot(tp2, scaled_pGcrA2, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('GcrA')
% xlabel('Time (min)')
% ylabel('Scaled GcrA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/GcrA_protein_full.eps','Resolution',300)
% end

% SciP protein plot
%scaling data points:
scaled_pSciP2 = (pSciP2 - min(pSciP2))/(max(pSciP2)-min(pSciP2))*(max(yout(a:b, SciP))-min(yout(a:b, SciP)))+min(yout(a:b, SciP));

figure()
 
line(tout(a:b), yout(a:b, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp2, scaled_pSciP2, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation', 'Empiracle')
xlim([1200 1350])
title('SciP')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/SciP_protein.eps','Resolution',300)
end

%SciP protein full time plot
% figure()
%  
% line(tout, yout(:, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tp2, scaled_pSciP2, 'ro', 'MarkerFaceColor', 'r')
% xlim([0 1500])
% title('SciP')
% xlabel('Time (min)')
% ylabel('Scaled SciP')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/SciP_protein_full.eps','Resolution',300)
% end

% CtrA protein plot
%scaling data:
scaled_pCtrA3 = (pCtrA3 - min(pCtrA3))/(max(pCtrA3)-min(pCtrA3))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA));
scaled_pCtrA5 = (pCtrA5 - min(pCtrA5))/(max(pCtrA5)-min(pCtrA5))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA));

figure()
 
line(tout(a:b), yout(a:b, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp3, scaled_pCtrA3, 'ro', 'MarkerFaceColor', 'r')
plot(tp5, scaled_pCtrA5, 'b^', 'MarkerFaceColor', 'b')

legend('Simulation','Empiracle Source 1', 'Empiracle Source 2', 'location', 'northwest')
xlim([1200 1350])
title('CtrA')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CtrA_protein.eps','Resolution',300)
end

% CtrA protein full time plot
% figure()
%  
% line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tp3, scaled_pCtrA3, 'ro', 'MarkerFaceColor', 'r')
% plot(tp5, scaled_pCtrA5, 'b^', 'MarkerFaceColor', 'b')
% xlim([0 1500])
% title('CtrA')
% xlabel('Time (min)')
% ylabel('Scaled CtrA')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/CtrA_protein_full.eps','Resolution',300)
% end
%%

% Complex3 protease plot
%figure(11)
%line(tout, yout(:, RcdA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
%hold on;
% plot(tp, pCtrA1, 'r*-')
%legend('Simulation')
%title('CPLX3')

% scatter(t, dCtrA, 'k*');
% scatter(t, dDnaA, 'm*');
% scatter(t, dGcrA, 'b*');
% scatter(t, dCcrM, 'y*');
% scatter(t, dSciP, 'g*');
% legend('Simulation','Experiment')
% title('mRNA of CtrA')
% axis([0 450 0 8])

%% Plotting DNA synthesis and Methylation variables

figure()
 
box on;
title('Probability of Hemimethylated States')
line(tout(a:b), yout(a:b, hcori), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p2 = line(tout(a:b), yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout(a:b), yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
h = legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}', 'Location', 'North');
xlabel('Time (min)')
ylabel('Probability')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/MethylationVars.eps','Resolution',300)
end

figure()
title('DNA synthesis')
box on;
p2 = line(tout(a:b), yout(a:b, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout(a:b), yout(a:b, DNA), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p4 = line(tout(a:b), yout(a:b, Count), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
% p5 = line(tout, yout(:, Ini), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
h = legend( 'Elongation',  'DNA', 'Chromosome');
xlabel('Time (min)')
ylabel('Count')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DNASynth.eps','Resolution',300)
end


%%  Protein degredation (model II) figures 
% plot Complex 1
figure();
set(gcf,'Position',[100 100 500 500])
hold on; 
box on;
plot(tout(a:b), yout(a:b, CPLX1), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([1200 1350]);
xlabel('Time (min)')
ylabel('Normalized Concentration')
title('Complex 1')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex1.eps','Resolution',300)
end

% % plot Complex 1 full time
% figure();
% set(gcf,'Position',[100 100 500 500])
% hold on; 
% box on;
% plot(tout, yout(:, CPLX1), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% xlim([0 1500]);
% xlabel('Time (min)')
% ylabel('Scaled Complex1')
% title('Complex 1')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/Complex1_full.eps','Resolution',300)
% end

%plot clpdrp
%scaled data:
scaled_pCpdRP = (pCpdRP - min(pCpdRP))/(max(pCpdRP)-min(pCpdRP))*(max(yout(a:b, CpdRP))-min(yout(a:b, CpdRP)))+min(yout(a:b, CpdRP));

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, CpdRP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpCpdRP, scaled_pCpdRP, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Empiracle')
xlim([1200 1350])
title('CpdR~P')
xlabel('Time (min)')
ylabel('Normalized Concentration')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP_protein.eps','Resolution',300)
end

% % plot CpdRP full time
% figure();
% set(gcf,'Position',[100 100 500 500])
% hold on;
% box on;
% plot(tout, yout(:, CpdRP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% plot(tpCpdRP, scaled_pCpdRP, 'ro', 'MarkerFaceColor', 'r');
% xlim([0 1500])
% xlabel('Time (min)')
% ylabel('Scaled CpdRP')
% title('CpdRP')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/CpdRP_protein_full.eps','Resolution',300)
% end

% plot total CpdR + CpdRP
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b, CpdRP) + yout(a:b, CpdR), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized Concentration')
title('Total CpdR (CpdR + CpdRP)')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP.eps','Resolution',300)
end

% % plot total CpdR + CpdRP full time
% figure();
% set(gcf,'Position',[100 100 500 500])
% hold on;
% box on;
% plot(tout, yout(:, CpdRP) + yout(:, CpdR), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% xlim([0 1500])
% xlabel('Time (min)')
% ylabel('Scaled Total CpdRP')
% title('Total CpdR (CpdR + CpdRP)')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/CpdRP_full.eps','Resolution',300)
% end

% plot Complex 2
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b, CPLX2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Concentration')
title('Complex 2')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex2.eps','Resolution',300)
end


% % plot Complex 2 full time
% figure();
% set(gcf,'Position',[100 100 500 500])
% hold on;
% box on;
% plot(tout, yout(:, CPLX2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% xlabel('Time (min)')
% ylabel('Scaled Complex2')
% title('Complex 2')
% xlim([0 1500])
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/Complex2_full.eps','Resolution',300)
% end

% plot Complex 3
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b,CPLX3), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Concentration')
title('Complex 3')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex3.eps','Resolution',300)
end


% plot Complex 3 full time
% figure();
% set(gcf,'Position',[100 100 500 500])
% hold on;
% box on;
% plot(tout, yout(:,CPLX3), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% xlabel('Time (min)')
% ylabel('Scaled Complex3')
% title('Complex 3')
% xlim([0 1500])
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/Complex3_full.eps','Resolution',300)
% end

% plot CpdR
time=[0,20,40,60,80,100,120] + 1200;
figure()
set(gcf,'Position',[100 100 500 500])
title('CpdR')
hold on;
box on;
CPDR = yout(:, CpdR);                                   
CPDR = CPDR(a:b);                  % gathering relevant simulated CpdR data
% CPDR=CPDR/max(CPDR);
%DIF = max(CPDR) - min(CPDR);
%CPDR = (CPDR-min(CPDR))/DIF;

plot(tout(a:b), CPDR, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')      % plotting simulated CpdR
% cpdr=[21001.196,20613.125,13581.933,10400.397,10563.811,13216.569,20216.276];
cpdr=[27532.359,40939.622,20027.844,10400.397,10563.811,13216.569,20216.276];
cpdr = (cpdr - min(cpdr))/(max(cpdr)-min(cpdr))*(max(CPDR)-min(CPDR))+min(CPDR);

% cpdr=cpdr/max(cpdr);

plot(time+10,cpdr,'ro','MarkerFaceColor','r')  % plotting experimental cpdr points

xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized Concentration')
legend('Simulation','Empiracle')
hold on;
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdR.eps','Resolution',300)
end

% % plot CpdR full time
% time=[0,20,40,60,80,100,120];
% figure()
% set(gcf,'Position',[100 100 500 500])
% title('CpdR')
% hold on;
% box on;
% CPDR = yout(:, CpdR);    
% 
% DIF = max(CPDR(1200:1350)) - min(CPDR(1200:1350));
% CPDR = (CPDR-min(CPDR))/DIF;
% plot(tout, CPDR, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')      % plotting simulated CpdR
% cpdr = (cpdr - min(cpdr))/(max(cpdr)-min(cpdr))*(max(CPDR(a:b))-min(CPDR(a:b)))+min(CPDR(a:b));
% 
% scatter(time+1210, cpdr, 'ro','MarkerFaceColor','r')  % plotting experimental cpdr points
% xlabel('Time (min)')
% ylabel('Scaled CpdR')
% hold on;
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/CpdR_full.eps','Resolution',300)
% end

% plot RcdA
figure()
set(gcf,'Position',[100 100 500 500])
title('RcdA')
hold on;
box on;
RCDA = yout(:, RcdA);
RCDA = RCDA(a:b);                    % gathering relevant simulated CpdR data
plot(tout(a:b), RCDA, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
xlim([1200 1350])

rcda = (rcda - min(rcda))/(max(rcda)-min(rcda))*(max(RCDA)-min(RCDA))+min(RCDA);

plot(tpRcdA, rcda, 'ro', 'MarkerFaceColor', 'r')  % plotting experimental rcda points
xlabel('Time (min)')
ylabel('Normalized Concentration')
legend('Simulation', 'Empiracle')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/RcdA.eps','Resolution',300)
end

%% plotting cdG and DivK~P
% cdG plot
% scaling data:
scaled_pcdG = (pcdG - min(pcdG))/(max(pcdG)-min(cdG(a:b)))*(max(cdG(a:b))-min(cdG(a:b)))+min(cdG(a:b));

figure()
 
line(tout(a:b) , cdG(a:b), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpcdG, scaled_pcdG, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
title('cdG')
xlabel('Time (min)')
ylabel('Normalized Concentration')
legend('Simulation', 'Empiracle')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/cdG.eps','Resolution',300)
end

% plotting cdG full time:
% figure()
%  
% line(tout , cdG, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tpcdG, scaled_pcdG, 'ro', 'MarkerFaceColor', 'r')
% xlim([0 1500])
% title('cdG')
% xlabel('Time (min)')
% ylabel('Normalized Concentration')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/cdG_full.eps','Resolution',300)
% end

% DivKP plot
% scaling data:
scaled_pDivKP = (pDivKP - min(pDivKP))/(max(pDivKP)-min(DivKP(a:b)))*(max(DivKP(a:b))-min(DivKP(a:b)))+min(DivKP(a:b));

figure()
 
line(tout(a:b) , DivKP(a:b), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpDivKP, scaled_pDivKP, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
title('DivK~P')
xlabel('Time (min)')
ylabel('Normalized Concentration')
legend('Simulation', 'Empiracle', 'location', 'northwest')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DivKP.eps','Resolution',300)
end

% plotting DivK~P full time:
% 
% figure()
%  
% line(tout , DivKP, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% hold on;
% box on;
% plot(tpDivKP, scaled_pDivKP, 'ro', 'MarkerFaceColor', 'r')
% xlim([0 1500])
% title('DivK~P')
% xlabel('Time (min)')
% ylabel('Scaled DivK~P')
% legend('Simulated DivK~P', 'Scaled Experimental Data')
% f = gcf;
% if save_figs == 1
%     exportgraphics(f,'./resources/generated_plots/DivKP_full.eps','Resolution',300)
% end
%% Old plotting code

% subplot(3,1,1);
% % cpdrp=[9118.983,4425.406,3905.163,1912.92,5768.719,8952.497,6521.347];
% cpdrp=[9118.983;4425.406;3905.163;1912.92;5768.719;8952.497;22145.296];
% % cpdrp=cpdrp/max(cpdrp);
% dif=max(cpdrp)-min(cpdrp);
% cpdrp=(cpdrp-min(cpdrp))/dif;
% 
% scatter(time+20,cpdrp,'ro','MarkerFaceColor','r')
% hold on;
% CPDRP=yout(:,3);
% % plot(tout,CPDRP);
% CPDRP=CPDRP(a:b);
% DIF=max(CPDRP)-min(CPDRP);
% CPDRP=(CPDRP-min(CPDRP))/DIF;
% plot(tout(a:b)-450,CPDRP)%plot the forth cell cycle 450-600min
% xlabel('Time/min')
% ylabel('CpdRP')
% legend('experimental data','simulated CpdRP')
% hold on;

% figure(55)
% line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% axis([0 420 0 0.1])

% figure(66)
% % subplot(2,3,1)
% figure(61)
% avg = mean(yout(:, CcrM))*ones(1,30);
% space1 = zeros(1,90);
% space2 = zeros(1,30);
% et = 1:1:150;
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CcrM)>mean(yout(:, CcrM)))=mean(yout(:, CcrM));
% scatter(tout, ybar, '.b', 'LineWidth', 4);
% plot(tout, yout(:,CcrM), 'k')
% title('CcrM')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 1])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% 
% % subplot(2,3,2)
% figure(62)
% avg = mean(yout(:, DnaA))*ones(1,50);
% space1 = zeros(1,10);
% space2 = zeros(1,90);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, DnaA)>mean(yout(:, DnaA)))=mean(yout(:, DnaA));
% scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,DnaA), 'k')
% ('DnaA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 3])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,3)
% figure(63)
% avg = mean(yout(:, GcrA))*ones(1,55);
% space1 = zeros(1,45);
% space2 = zeros(1,50);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, GcrA)>mean(yout(:, GcrA)))=mean(yout(:, GcrA));
% scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,GcrA), 'k')
% title('GcrA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 6])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,4)
% figure(64)
% avg = mean(yout(:, SciP))*ones(1,15);
% avg2 = mean(yout(:, SciP))*ones(1,60);
% space1 = zeros(1,85);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, SciP)>mean(yout(:, SciP)))=mean(yout(:, SciP));
% scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,SciP), 'k')
% title('SciP')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 7])
% grid on 
% set(gca,'xtick',[0:30:150])
% 
% % subplot(2,3,5)
% figure(65)
% avg = mean(yout(:, CtrA))*ones(1,15);
% avg2 = mean(yout(:, CtrA))*ones(1,60);
% space1 = zeros(1,75);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'r*', 'LineWidth', 4);
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CtrA)>mean(yout(:, CtrA)))=mean(yout(:, CtrA));
% scatter(tout, ybar, 'b.', 'LineWidth', 4);
% plot(tout, yout(:,CtrA), 'k')
% title('CtrA')
% legend('Experiment','Simulation','Curve')
% axis([0 150 0 8])

% grid on 
% set(gca,'xtick',[0:30:150])


% h = legend('CtrA', 'DnaA',  'GcrA', 'CcrM', 'Location', 'North');
% axis([0 500 0.2 4])

% set(gca,'YTick',[1:4]) 
% lab = ['CcrM';'GcrA'; 'DnaA';'CtrA'];
% set(gca,'yticklabel',lab);
% title('Temporal control results from Shenghua model')
%%


%end of plotting
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

% figure(66)
% subplot(5,1,1)
% % figure(61)
% avg = ones(1,30);
% space1 = zeros(1,90);
% space2 = zeros(1,30);
% et = 1:1:150;
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CcrM)>mean(yout(:, CcrM)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('CcrM')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% 
% subplot(5,1,2)
% % figure(62)
% avg = ones(1,50);
% space1 = zeros(1,10);
% space2 = zeros(1,90);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, DnaA)>mean(yout(:, DnaA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('DnaA')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,3)
% % figure(63)
% avg = ones(1,55);
% space1 = zeros(1,45);
% space2 = zeros(1,50);
% e1 = [space1,  avg,  space2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, GcrA)>mean(yout(:, GcrA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('GcrA')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,4)
% % figure(64)
% avg = ones(1,15);
% avg2 = ones(1,60);
% space1 = zeros(1,75);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, SciP)>mean(yout(:, SciP)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('SciP')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})
% 
% subplot(5,1,5)
% % figure(65)
% avg = ones(1,15);
% avg2 = ones(1,90);
% space1 = zeros(1,45);
% e1 = [avg, space1, avg2];
% scatter(et, e1, 'or', 'LineWidth', 1, 'MarkerFaceColor', 'r' );
% hold on;
% ybar = zeros(size(tout));
% ybar(yout(:, CtrA)>mean(yout(:, CtrA)))=2;
% scatter(tout, ybar, 'sb', 'LineWidth', 2, 'MarkerFaceColor', 'b' );
% title('CtrA')
% % legend('Experiment','Simulation','Curve')
% axis([0 150 0.5 2.5])
% grid on 
% set(gca,'xtick',[0:30:150])
% set(gca,'ytick',[1,2])
% set(gca,'yticklabel',{'Experiment','Simulation'})



