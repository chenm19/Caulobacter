%% Ploting the time courses of model variables
close all; clf;

load('output.mat')

save_figs = 0;  %switch for saving: 1 = save; 0 = no save

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Five transcriptional regulators
% figure(11)
% title('Four transcriptional regulators from Shenghua model')
% p1 = line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
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

% old protein data
tp = [0, 20, 40, 60, 80, 100, 120, 140] + 1200;
tpCcrM = [0, 20, 40, 60, 80, 100, 120] + 1200;
tp2_old = [9, 27, 44, 62, 80, 98, 117, 134, 152] + 1200;

pDnaA = 2*[0.4567	0.9857	0.9147	0.3154	0.128	0.0793	0.4203	0.6582];
pGcrA = 6*[0	0.5021	1	0.884	1	0.554	0.662	0.7388];
pCtrA1_old = 7*[0.7791	0	0	0	0.4682	1	0.8591	0.6672];
pCcrM = [0.979012341	0.245634747	0.105451343	0.054014407	0.031971569	0.078998499	1];
pSciP = 0.7*[7.757911401	7.931725156	4.139518022	1	0.446070124	0.429447455	0.761125117	3.856508369	8.881823903];
pCtrA2_old = [0.8158	0.614	0.1881	0.0323	0.3822	0.7033	0.8924	1	0.782];

%% new protein data:

% publication 1:
tp1 = [0 20 40 60 80 100 120 140]
pDnaA1 = [4784.305 10537.962 9511.012 3134.012 1447.184 1014.355 3693.083 6339.548]
pGcrA1 = [0 8632.062 16136.355 14185.012 14835.891 8353.184 9923.598 10884.719]
pCtrA1 = [6648.012 0 0 0 4042.598 8955.376 7497.598 6225.962]

%publication 2:
tp2 = [0 20 40 60 80 100 120 140 160]
pSciP2 = [9239.255 9518.134 4585.305 610.042 161.657 388.435 4424.062 10748.426 12927.062]
pGcrA2 = [0 444.87 4548.648 6148.669 5947.891 4592.062 5067.406 5809.548 6002.77]
pCtrA2 = [7177.962 9270.255 1045.406 0 315.506 5709.77 8172.598 7724.355 8002.941]

%publication 3:
tp3 = [0 10 20 30 40 50 60 70 80 90 100]
pCcrM3 = [1202.5475 1202.5475 1861.891 2587.134 3223.598 4278.033 3304.477 6183.619 6666.841 5629.205]
pCtrA3 = [4366.669 544.698 678.284 1565.962 3814.74 6857.891 7663.376 13218.74 12383.598 10310.083]

%publication 4:
tp4 = tp3
pSciP4 = [8176.062 3082.648 894.456 248.12 248.12 248.12 248.12 855.77 1727.527 2709.941]
pCtrA4 = [3465.406 1278.77 928.527 590.991 1960.355 3866.477 4892.77 5509.062 4213.941 4028.77]


%%
% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350));

%% Plotting concentration of 5 main regulators mRNA against time

% ccrM mRNA plot
figure(1)
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCcrM, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
legend('Simulated ccrm','Experimental data', 'location', 'northeastoutside')
title('\it{ccrM} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized ccrM')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ccrM.png','Resolution',300)
end

% ccrM mRNA full plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCcrM, 'ro', 'MarkerFaceColor', 'r')
title('\it{ccrM}')
xlabel('Time (min)')
ylabel('Normalized ccrM')
f = gcf;
xlim([0 1500])
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ccrM_full.png','Resolution',300)
end

% dnaA mRNA plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dDnaA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulated dnaA','Experimental data', 'location', 'northeastoutside')
title('\it{dnaA} (8th cell cycle)')
xlabel('Time (min)')
xlim([1200 1350])
ylabel('Normalized dnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/dnaA.png','Resolution',300)
end

% plot dnaA full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dDnaA, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('\it{dnaA}')
xlabel('Time (min)')
xlim([0 1500])
ylabel('Normalized dnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/dnaA_full.png','Resolution',300)
end


% gcra mRNA plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dGcrA, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
legend('Simulated gcra','Experimental data', 'location', 'northeastoutside')
title('\it{gcrA} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized gcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/gcrA.png','Resolution',300)
end

% gcra mRNA full time subplot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dGcrA, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('\it{gcrA}')
xlabel('Time (min)')
ylabel('Normalized gcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/gcrA_full.png','Resolution',300)
end

% sciP mRNA plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dSciP, 'ro', 'MarkerFaceColor', 'r')
legend('Simulated sciP','Experimental data', 'location', 'northeastoutside')
xlim([1200 1350])
title('\it{sciP} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized sciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/sciP.png','Resolution',300)
end

% sciP mRNA full time subplot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dSciP, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('\it{sciP}')
xlabel('Time (min)')
ylabel('Normalized sciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/sciP_full.png','Resolution',300)
end

% ctra mRNA plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCtrA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulated ctrA','Experimental data', 'location', 'northeastoutside')
xlim([1200 1350])
title('\it{ctrA} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized crrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ctrA.png','Resolution',300)
end

% ctra mRNA full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCtrA, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('\it{ctrA}')
xlabel('Time (min)')
ylabel('Normalized crrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ctrA_full.png','Resolution',300)
end

%% Plotting concentration of 5 main regulator proteins (and protease complex) against time
% Ccrm protein plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpCcrM, pCcrM, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Experiment', 'location', 'northeastoutside')
xlim([1200 1350])
title('CcrM (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized CcrM')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Ccrm_protein.png','Resolution',300)
end

% Ccrm protein full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpCcrM, pCcrM, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('CcrM')
xlabel('Time (min)')
ylabel('Normalized CcrM')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Ccrm_protein_full.png','Resolution',300)
end

% DnaA protein plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pDnaA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Experiment', 'location', 'northeastoutside')
xlim([1200 1350])
title('DnaA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized DnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DnaA_protein.png','Resolution',300)
end

% DnaA protein full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pDnaA, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('DnaA')
xlabel('Time (min)')
ylabel('Normalized DnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DnaA_protein_full.png','Resolution',300)
end

% GcrA protein plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pGcrA, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Experiment', 'location', 'northeastoutside')
xlim([1200 1350])
title('GcrA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized GcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/GcrA_protein.png','Resolution',300)
end

% GcrA protein full time subplot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pGcrA, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('GcrA')
xlabel('Time (min)')
ylabel('Normalized GcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/GcrA_protein_full.png','Resolution',300)
end

% SciP protein plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp2, pSciP, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Experiment', 'location', 'northeastoutside')
xlim([1200 1350])
title('SciP (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized SciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/SciP_protein.png','Resolution',300)
end

% SciP protein full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp2, pSciP, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('SciP')
xlabel('Time (min)')
ylabel('Normalized SciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/SciP_protein_full.png','Resolution',300)
end

% CtrA protein plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pCtrA1, 'ro', 'MarkerFaceColor', 'r')
legend('Simulation','Experiment', 'location', 'northeastoutside')
xlim([1200 1350])
title('CtrA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Normalized CtrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CtrA_protein.png','Resolution',300)
end

% CtrA protein full time plot
figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout, yout(:, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp, pCtrA1, 'ro', 'MarkerFaceColor', 'r')
xlim([0 1500])
title('CtrA')
xlabel('Time (min)')
ylabel('Normalized CtrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CtrA_protein_full.png','Resolution',300)
end
%%

% Complex3 protease plot
%figure(11)
%p1 = line(tout, yout(:, RcdA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
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
set(gcf,'Position',[100 100 500 500])
box on;
title('Probability of Hemimethylated States')
p1 = line(tout(a:b), yout(a:b, hcori), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p2 = line(tout(a:b), yout(a:b, hCcrM), 'Color', 'm', 'LineWidth', 2, 'Linestyle', '-');
p3 = line(tout(a:b), yout(a:b, hCtrA), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
h = legend('h_{Cori}', 'h_{ccrM}',  'h_{ctrA}', 'Location', 'North');
xlabel('Time (min)')
ylabel('Probability')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/MethylationVars.png','Resolution',300)
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
    exportgraphics(f,'./resources/generated_plots/DNASynth.png','Resolution',300)
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
ylabel('Normalized Complex1')
title('Complex 1 (8th cell cycle)')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex1.png','Resolution',300)
end

% plot Complex 1 full time
figure();
set(gcf,'Position',[100 100 500 500])
hold on; 
box on;
plot(tout, yout(:, CPLX1), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([0 1500]);
xlabel('Time (min)')
ylabel('Normalized Complex1')
title('Complex 1')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex1_full.png','Resolution',300)
end

% plot CpdRP
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b, CpdRP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized CpdRP')
title('CpdRP (8th cell cycle)')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP.png','Resolution',300)
end

% plot CpdRP full time
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout, yout(:, CpdRP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([0 1500])
xlabel('Time (min)')
ylabel('Normalized CpdRP')
title('CpdRP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP_full.png','Resolution',300)
end

% plot total CpdR + CpdRP
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b, CpdRP) + yout(a:b, CpdR), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized CpdRP')
title('Total CpdR (CpdR + CpdRP) (8th cell cycle)')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP.png','Resolution',300)
end

% plot total CpdR + CpdRP full time
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout, yout(:, CpdRP) + yout(:, CpdR), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlim([0 1500])
xlabel('Time (min)')
ylabel('Normalized CpdRP')
title('Total CpdR (CpdR + CpdRP)')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdRP_full.png','Resolution',300)
end

% plot Complex 2
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b, CPLX2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Complex2')
title('Complex 2 (8th cell cycle)')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex2.png','Resolution',300)
end

% plot Complex 2 full time
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout, yout(:, CPLX2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Complex2')
title('Complex 2')
xlim([0 1500])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex2_full.png','Resolution',300)
end

% plot Complex 3
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout(a:b), yout(a:b,CPLX3), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Complex3')
title('Complex 3 (8th cell cycle)')
xlim([1200 1350])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex3.png','Resolution',300)
end


% plot Complex 3 full time
figure();
set(gcf,'Position',[100 100 500 500])
hold on;
box on;
plot(tout, yout(:,CPLX3), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
xlabel('Time (min)')
ylabel('Normalized Complex3')
title('Complex 3')
xlim([0 1500])
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Complex3_full.png','Resolution',300)
end

% plot CpdR
time=[0,20,40,60,80,100,120] + 1200;
figure()
set(gcf,'Position',[100 100 500 500])
title('CpdR (8th cell cycle)')
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

scatter(time+10,cpdr,'ro','MarkerFaceColor','r')  % plotting experimental cpdr points
xlim([1200 1350])
xlabel('Time (min)')
ylabel('Normalized CpdR')
legend('experimental data', 'simulated CpdR', 'location', 'northeastoutside')
hold on;
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdR.png','Resolution',300)
end

% plot CpdR full time
time=[0,20,40,60,80,100,120];
figure()
set(gcf,'Position',[100 100 500 500])
title('CpdR')
hold on;
box on;
CPDR = yout(:, CpdR);    

DIF = max(CPDR(1200:1350)) - min(CPDR(1200:1350));
CPDR = (CPDR-min(CPDR))/DIF;
plot(tout, CPDR, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')      % plotting simulated CpdR
cpdr = (cpdr - min(cpdr))/(max(cpdr)-min(cpdr))*(max(CPDR(a:b))-min(CPDR(a:b)))+min(CPDR(a:b));

scatter(time+1210, cpdr, 'ro','MarkerFaceColor','r')  % plotting experimental cpdr points
xlabel('Time (min)')
ylabel('Normalized CpdR')
hold on;
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CpdR_full.png','Resolution',300)
end

% plot RcdA
figure()
set(gcf,'Position',[100 100 500 500])
title('RcdA (8th cell cycle)')
hold on;
box on;
RCDA = yout(:, RcdA);
RCDA = RCDA(a:b);                    % gathering relevant simulated CpdR data
% RCDA=RCDA/max(RCDA);
%DIF = max(RCDA) - min(RCDA);
%RCDA = (RCDA-min(RCDA))/DIF;         % plotting simulated RcdA data

plot(tout(a:b), RCDA, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
xlim([1200 1350])
time2=[0,9.4,26.7,44.8,62.1176,79.8118,97.5059,115.2,133.271,151.341];
rcda=[0.0000136,0.093,0.984,1.00163,0.528269,0.517348,0.517987,0.431921,0.738932,0.271377];
rcda = (rcda - min(rcda))/(max(rcda)-min(rcda))*(max(RCDA)-min(RCDA))+min(RCDA);
scatter(time2 + 1200, rcda, 'ro', 'MarkerFaceColor', 'r')  % plotting experimental rcda points
xlabel('Time (min)')
ylabel('Normalized RcdA')
legend('simulated RcdA', 'experimental data', 'location', 'northeastoutside')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/RcdA.png','Resolution',300)
end

% plot RcdA full time
figure()
set(gcf,'Position',[100 100 500 500])
title('RcdA')
hold on;
box on;
RCDA = yout(:, RcdA);
plot(tout, RCDA, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-')
xlim([0 1500])

scatter(time2 + 1200, rcda, 'ro', 'MarkerFaceColor', 'r')  % plotting experimental rcda points
xlabel('Time (min)')
ylabel('Normalized RcdA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/RcdA_full.png','Resolution',300)
end

%% plotting cdG and DivK~P
T = 150;        % period of Caulobacter
t_d = rem(tout,T); % return remainder after division t/T

% cdG plot
% data
cdG= 0.3233 * sin(pi/75 * t_d + 0.3469) + 0.3363;
cdG(cdG<0) = 0;
time = [0    20    40    60    80   100   120]+1200;
cdG_levels = [0.2292 1.0000 0.5000 0.2833 0.2250 0.1625 0];

figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b) , cdG(a:b), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(time, cdG_levels, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
title('cdG (8th Cell Cycle)')
xlabel('Time (min)')
ylabel('Normalized cdG')
legend('Simulated cdG', 'Scaled Experimental Data', 'location', 'northeastoutside')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/cdG.png','Resolution',300)
end

% DivKP plot
% data
DivKP = 0.004218 * t_d + 0.5779;
DivKP_levels = [0.6311 0.7972 0.8103 0.916 1]
time = [20 40 60 80 100] + 1200;

figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b) , DivKP(a:b), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(time, DivKP_levels, 'ro', 'MarkerFaceColor', 'r')
xlim([1200 1350])
title('DivK~P (8th Cell Cycle)')
xlabel('Time (min)')
ylabel('Normalized DivK~P')
legend('Simulated DivK~P', 'Scaled Experimental Data', 'location', 'northeastoutside')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DivKP.png','Resolution',300)
end
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
% p1 = line(tout, yout(:, Elong), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
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
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
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
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
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
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
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
% p1 = scatter(tout, ybar, 'b.', 'LineWidth', 4);
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



