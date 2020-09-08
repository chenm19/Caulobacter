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

% isolating the indices of the 8th cell cycle:
[~, a]=min(abs(tout(:)-1200));  % a and b are indeces of beginning and end of cell cycle 
[~, b]=min(abs(tout(:)-1350)); 

%% bar chart test
figure(1)
box on;
X = categorical({'DnaA','GcrA','CtrA','CcrM', 'SciP'});
max_DnaA = max(yout(a:b, DnaA));
max_GcrA = max(yout(a:b, GcrA));
max_CtrA = max(yout(a:b, CtrA));
max_CcrM = max(yout(a:b, CcrM));
max_SciP = max(yout(a:b, SciP));
Y = [max_DnaA, max_GcrA, max_CtrA, max_CcrM, max_SciP];
bar(X,Y,'FaceColor',[0 0.4470 0.7410],'EdgeColor',[1 1 1])
title('Relative Maximum Protein Concentrations')
if save_figs == 1
     exportgraphics(f,'./resources/generated_plots/bar_chart.eps','Resolution',300)
end
%% line chart

figure(2)
subplot(5,1,1);
avg = mean(yout(a:b, CcrM));
simulation_abv_avg = (yout(a:b, CcrM)) > avg;

scatter(tout(a:b)-1200, simulation_abv_avg, 'MarkerFaceColor', 'r');
hold on;
ybar = zeros(150);
ybar(yout(:, CcrM) > mean(yout(:, CcrM))) = mean(yout(:, CcrM));
%scatter(tout(a:b), ybar, '^b', 'MarkerFaceColor', 'b');
%plot(tout, yout(:,CcrM), 'k')
title('CcrM')
legend('Experiment','Simulation','location', 'westoutside')
legend('boxoff')
axis([0 150 0 1])
grid on 
set(gca,'xtick',[0:30:150])


subplot(5, 1, 2);
avg = mean(yout(a:b, DnaA));
simulation_abv_avg = (yout(a:b, DnaA)) > avg;

scatter(tout(a:b)-1200, simulation_abv_avg, 'MarkerFaceColor', 'r');
hold on;
ybar = zeros(150);
ybar(yout(:, DnaA) > mean(yout(:, DnaA))) = mean(yout(:, DnaA));
%scatter(tout, ybar, '^b', 'MarkerFaceColor', 'b');
%plot(tout, yout(:,CcrM), 'k')
title('DnaA')
legend('Experiment','Simulation','location', 'westoutside')
legend('boxoff')
axis([0 150 0 1])
grid on 
set(gca,'xtick',[0:30:150])
 
subplot(5,1,3)

avg = mean(yout(a:b, GcrA));
simulation_abv_avg = (yout(a:b, GcrA)) > avg;

scatter(tout(a:b)-1200, simulation_abv_avg, 'MarkerFaceColor', 'r');
hold on;
ybar = zeros(150);
ybar(yout(:, GcrA) > mean(yout(:, GcrA))) = mean(yout(:, GcrA));
%scatter(tout, ybar, '^b', 'MarkerFaceColor', 'b');
title('GcrA')
legend('Experiment','Simulation','location', 'westoutside')
legend('boxoff')
axis([0 150 0 1])
grid on 
set(gca,'xtick',[0:30:150])

subplot(5,1,4)
avg = mean(yout(a:b, SciP));
simulation_abv_avg = (yout(a:b, SciP)) > avg;

scatter(tout(a:b)-1200, simulation_abv_avg, 'MarkerFaceColor', 'r');
hold on;
ybar = zeros(150);
ybar(yout(:, SciP) > mean(yout(:, SciP))) = mean(yout(:, SciP));
%scatter(tout, ybar, '^b', 'MarkerFaceColor', 'b');
%plot(tout, yout(:,CcrM), 'k')
title('SciP')
legend('Experiment','Simulation','location', 'westoutside')
legend('boxoff')
axis([0 150 0 1])
grid on 
set(gca,'xtick',[0:30:150])

subplot(5,1,5)
avg = mean(yout(a:b, CtrA));
simulation_abv_avg = (yout(a:b, CtrA)) > avg;

scatter(tout(a:b)-1200, simulation_abv_avg, 'MarkerFaceColor', 'r');
hold on;
ybar = zeros(150);
ybar(yout(:, CtrA) > mean(yout(:, CtrA))) = mean(yout(:, CtrA));
%scatter(tout, ybar, '^b', 'MarkerFaceColor', 'b');
%plot(tout, yout(:,CcrM), 'k')
title('CtrA')
legend('Experiment','Simulation','location', 'westoutside')
legend('boxoff')
axis([0 150 0 1])
grid on 
set(gca,'xtick',[0:30:150])