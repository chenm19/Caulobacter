%% Ploting the time courses of model variables
close all; clf;

load('output.mat')

save_figs = 0;  %switch for saving: 1 = save; 0 = no save


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


%%
% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350));

%% Plotting concentration of 5 main regulators mRNA against time

% ccrM mRNA plot
% scaling data:
scaled_dCcrM = (dCcrM - min(dCcrM))/(max(dCcrM)-min(dCcrM))*(max(yout(a:b, mCcrM))-min(yout(a:b, mCcrM)))+min(yout(a:b, mCcrM));

figure(1)
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); %simulated levels
hold on;
box on;
plot(t, dCcrM, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dCcrM, 'b^', 'MarkerFaceColor', 'b')
xlim([1200 1350])
legend('Simulated ccrm','Experimental Data', 'location', 'northeastoutside')
title('\it{ccrM} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled ccrM')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ccrM.png','Resolution',300)
end

% ccrM mRNA full plot
figure()
set(gcf,'Position',[100 100 500 500])
line(tout, yout(:, mCcrM), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCcrM, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dCcrM, 'b^', 'MarkerFaceColor', 'b')
title('\it{ccrM}')
xlabel('Time (min)')
ylabel('Scaled ccrM')
f = gcf;
xlim([0 1500])
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ccrM_full.png','Resolution',300)
end

% dnaA mRNA plot
% scaled data:
scaled_dDnaA = (dDnaA - min(dDnaA))/(max(dDnaA)-min(dDnaA))*(max(yout(a:b, mDnaA))-min(yout(a:b, mDnaA)))+min(yout(a:b, mDnaA));

figure()
set(gcf,'Position',[100 100 500 500])
p1 = line(tout(a:b), yout(a:b, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-'); % simulated dnaA
hold on;
box on;
plot(t, dDnaA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dDnaA, 'b^', 'MarkerFaceColor', 'b')
legend('Simulated dnaA','Experimental Data', 'location', 'northeastoutside')
title('\it{dnaA} (8th cell cycle)')
xlabel('Time (min)')
xlim([1200 1350])
ylabel('Scaled dnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/dnaA.png','Resolution',300)
end

% plot dnaA full time plot
figure()
set(gcf,'Position',[100 100 500 500])
line(tout, yout(:, mDnaA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dDnaA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dDnaA, 'b^', 'MarkerFaceColor', 'b')
xlim([0 1500])
title('\it{dnaA}')
xlabel('Time (min)')
xlim([0 1500])
ylabel('Scaled dnaA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/dnaA_full.png','Resolution',300)
end


% gcra mRNA plot
% scaled data:
scaled_dGcrA = (dGcrA - min(dGcrA))/(max(dGcrA)-min(dGcrA))*(max(yout(a:b, mGcrA))-min(yout(a:b, mGcrA)))+min(yout(a:b, mGcrA));

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dGcrA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dGcrA, 'b^', 'MarkerFaceColor', 'b')
xlim([1200 1350])
legend('Simulated gcra','Experimental Data', 'location', 'northeastoutside')
title('\it{gcrA} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled gcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/gcrA.png','Resolution',300)
end

% gcra mRNA full time subplot
figure()
set(gcf,'Position',[100 100 500 500])
line(tout, yout(:, mGcrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dGcrA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dGcrA, 'b^', 'MarkerFaceColor', 'b')
xlim([0 1500])
title('\it{gcrA}')
xlabel('Time (min)')
ylabel('Scaled gcrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/gcrA_full.png','Resolution',300)
end

% sciP mRNA plot
% scaled data:
scaled_dSciP = (dSciP - min(dSciP))/(max(dSciP)-min(dSciP))*(max(yout(a:b, mSciP))-min(yout(a:b, mSciP)))+min(yout(a:b, mSciP));

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dSciP, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dSciP, 'b^', 'MarkerFaceColor', 'b')
legend('Simulated sciP','Experimental Data', 'location', 'northeastoutside')
xlim([1200 1350])
title('\it{sciP} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled sciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/sciP.png','Resolution',300)
end

% sciP mRNA full time subplot
figure()
set(gcf,'Position',[100 100 500 500])
line(tout, yout(:, mSciP), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dSciP, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dSciP, 'b^', 'MarkerFaceColor', 'b')
xlim([0 1500])
title('\it{sciP}')
xlabel('Time (min)')
ylabel('Scaled sciP')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/sciP_full.png','Resolution',300)
end

% ctrA mRNA plot
%scaled data:
scaled_dCtrA = (dCtrA - min(dCtrA))/(max(dCtrA)-min(dCtrA))*(max(yout(a:b, mCtrA))-min(yout(a:b, mCtrA)))+min(yout(a:b, mCtrA));

figure()
set(gcf,'Position',[100 100 500 500])
line(tout(a:b), yout(a:b, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCtrA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dCtrA, 'b^', 'MarkerFaceColor', 'b')
legend('Simulated ctrA','Experimental Data', 'location', 'northeastoutside')
xlim([1200 1350])
title('\it{ctrA} (8th cell cycle)')
xlabel('Time (min)')
ylabel('Scaled crrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ctrA.png','Resolution',300)
end

% ctra mRNA full time plot
figure()
set(gcf,'Position',[100 100 500 500])
line(tout, yout(:, mCtrA), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on;
box on;
plot(t, dCtrA, 'ro', 'MarkerFaceColor', 'r')
plot(t, scaled_dCtrA, 'b^', 'MarkerFaceColor', 'b')
xlim([0 1500])
title('\it{ctrA}')
xlabel('Time (min)')
ylabel('Scaled crrA')
f = gcf;
if save_figs == 1
    exportgraphics(f,'./resources/generated_plots/ctrA_full.png','Resolution',300)
end
