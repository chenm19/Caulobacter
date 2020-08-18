load('output.mat');
close all;
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
%% new protein data:

% publication 1:
tp1 = [0 20 40 60 80 100 120 140] + 1200;
pDnaA1 = [4784.305 10537.962 9511.012 3134.012 1447.184 1014.355 3693.083 6339.548];
pGcrA1 = [0 8632.062 16136.355 14185.012 14835.891 8353.184 9923.598 10884.719];
pCtrA1 = [6648.012 0 0 0 4042.598 8955.376 7497.598 6225.962];

%publication 2:
tp2 = [0 20 40 60 80 100 120 140 160] + 1200;
pSciP2 = [9239.255 9518.134 4585.305 610.042 161.657 388.435 4424.062 10748.426 12927.062];
pGcrA2 = [0 444.87 4548.648 6148.669 5947.891 4592.062 5067.406 5809.548 6002.77];
pCtrA2 = [7177.962 9270.255 1045.406 0 315.506 5709.77 8172.598 7724.355 8002.941];

%publication 3:
tp3 = [0 10 20 30 40 50 60 70 80 90] + 1200;
pCcrM3 = [1202.5475 1202.5475 1861.891 2587.134 3223.598 4278.033 3304.477 6183.619 6666.841 5629.205];
pCtrA3 = [4366.669 544.698 678.284 1565.962 3814.74 6857.891 7663.376 13218.74 12383.598 10310.083];

%publication 4:
tp4 = [0 10 20 30 40 50 60 70 80 90 100] + 1200;
pSciP4 = [8176.062 3082.648 894.456 248.12 248.12 248.12 248.12 855.77 1727.527 2709.941 2193.042];
pCtrA4 = [3465.406 1278.77 928.527 590.991 1960.355 3866.477 4892.77 5509.062 4213.941 4028.77 2851.042];

%% Plotting concentration of 5 main regulator proteins (and protease complex) against time
% Ccrm protein plot
PCCRM3 = (pCcrM3 - min(pCcrM3))/(max(pCcrM3)-min(pCcrM3))*(max(yout(a:b, CcrM))-min(yout(a:b, CcrM)))+min(yout(a:b, CcrM))
figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, CcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tpCcrM_old, pCcrM_old, 'ro-', 'MarkerFaceColor', 'r')
plot(tp3, PCCRM3, 'go-', 'MarkerFaceColor', 'g')
legend('Simulation','Original data','Pub3 data', 'location', 'northeastoutside')
xlim([1200 1350])
title('CcrM (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled CcrM')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/Ccrm_protein.png','Resolution',300)
end

% DnaA protein plot
PDNAA1 = (pDnaA1 - min(pDnaA1))/(max(pDnaA1)-min(pDnaA1))*(max(yout(a:b, DnaA))-min(yout(a:b, DnaA)))+min(yout(a:b, DnaA))

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, DnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp_old, pDnaA_old, 'ro-', 'MarkerFaceColor', 'r')
plot(tp1, PDNAA1, 'go-', 'MarkerFaceColor', 'g')
legend('Simulation','Original data', 'Pub 1 data', 'northeastoutside')
xlim([1200 1350])
title('DnaA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled DnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/DnaA_protein.png','Resolution',300)
end

% GcrA protein plot
PGCRA1 = (pGcrA1 - min(pGcrA1))/(max(pGcrA1)-min(pGcrA1))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA))
PGCRA2 = (pGcrA2 - min(pGcrA2))/(max(pGcrA2)-min(pGcrA2))*(max(yout(a:b, GcrA))-min(yout(a:b, GcrA)))+min(yout(a:b, GcrA))

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, GcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp_old, pGcrA_old, 'ro-', 'MarkerFaceColor', 'r')
plot(tp1, PGCRA1, 'go-', 'MarkerFaceColor', 'g')
plot(tp2, PGCRA2, 'bo-', 'MarkerFaceColor', 'b')

legend('Simulation','Original data', 'Pub 1','pub 2', 'northeastoutside')
xlim([1200 1350])
title('GcrA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled GcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/GcrA_protein.png','Resolution',300)
end


% SciP protein plot
PSCIP2 = (pSciP2 - min(pSciP2))/(max(pSciP2)-min(pSciP2))*(max(yout(a:b, SciP))-min(yout(a:b, SciP)))+min(yout(a:b, SciP))
PSCIP4 = (pSciP4 - min(pSciP4))/(max(pSciP4)-min(pSciP4))*(max(yout(a:b, SciP))-min(yout(a:b, SciP)))+min(yout(a:b, SciP))

figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, SciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp2_old, pSciP_old, 'ro-', 'MarkerFaceColor', 'r')
plot(tp2, PSCIP2, 'go-', 'MarkerFaceColor', 'g')
plot(tp4, PSCIP4, 'bo-', 'MarkerFaceColor', 'b')
legend('Simulation','Original Data', 'Pub 2', 'Pub 4', 'location', 'northeastoutside')
xlim([1200 1350])
title('SciP (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled SciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/SciP_protein.png','Resolution',300)
end


% CtrA protein plot
PCTRA1 = (pCtrA1 - min(pCtrA1))/(max(pCtrA1)-min(pCtrA1))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA))
PCTRA2 = (pCtrA2 - min(pCtrA2))/(max(pCtrA2)-min(pCtrA2))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA))
PCTRA3 = (pCtrA3 - min(pCtrA3))/(max(pCtrA3)-min(pCtrA3))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA))
PCTRA4 = (pCtrA4 - min(pCtrA4))/(max(pCtrA4)-min(pCtrA4))*(max(yout(a:b, CtrA))-min(yout(a:b, CtrA)))+min(yout(a:b, CtrA))

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, CtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(tp_old, pCtrA1_old, 'ro-', 'MarkerFaceColor', 'r')
plot(tp1, PCTRA1, 'go-', 'MarkerFaceColor', 'g')
plot(tp2, PCTRA2, 'bo-', 'MarkerFaceColor', 'b')
plot(tp3, PCTRA3, 'mo-', 'MarkerFaceColor', 'm')
plot(tp4, PCTRA4, 'co-', 'MarkerFaceColor', 'c')
legend('Simulation','Original','Pub 1','Pub 2','Pub 3','Pub 4', 'location', 'northeastoutside')
xlim([1200 1350])
title('CtrA (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled CtrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/CtrA_protein.png','Resolution',300)
end

