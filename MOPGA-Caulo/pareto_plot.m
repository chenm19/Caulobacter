%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% Script for ploting pareto front and calculate the sensitivity and goodness 
% of fit.
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
%% load the optimization result

load('./LoadMatrices/opt_result.mat')

set(groot, 'DefaultAxesFontSize', 23)
set(groot, 'DefaultFigurePosition', [200 50 600 600])
set(groot, 'DefaultFigureColormap',jet)
set(groot, 'DefaultFigureGraphicsSmoothing','off')
set(groot, 'DefaultAxesLineWidth', 1.2)

% set to 1 to save figures
save_figs = 0;

figure(1)
hold on;
plot(fval(:,1), fval(:,2), 'm+', 'MarkerSize',15, 'LineWidth', 3)
plot(fval(sensitiveID,1), fval(sensitiveID,2), 'ko', 'MarkerSize',15, 'LineWidth', 3)
plot(fval(9,1), fval(9,2), 'r*', 'MarkerSize',  15, 'LineWidth', 3)
xlabel('$f_1(\chi)$','interpreter','latex','FontSize',25)
ylabel('$f_2(\chi)$','interpreter','latex','FontSize',25)
title('Pareto front','FontSize',40)
box on
legend('Nondominated points','Sensitivity analysis','Best for model')
f = gcf;    pbaspect([1.3 1 1])

if save_figs == 1
    exportgraphics(f,'./Figures/Pareto_Plot.eps','Resolution',300)
end


%% sensitivity
mean(sensitiveF(:,1)./fval(sensitiveID,1))-1
mean(sensitiveF(:,2)./fval(sensitiveID,2))-1

%% goodness of fitting, rmse
rms(fval(sensitiveID,1))
rms(fval(sensitiveID,2))
