load('opt_result.mat')

% selID = [4 9 15 11];

figure(1)
hold on;
plot(fval(:,1), fval(:,2), 'm+', 'MarkerSize',10, 'LineWidth', 2)
plot(fval(sensitiveID,1), fval(sensitiveID,2), 'ko', 'MarkerSize',15, 'LineWidth', 3)
plot(fval(9,1), fval(9,2), 'r*', 'MarkerSize',10, 'LineWidth', 2)
xlabel('$f_1(\chi)$','interpreter','latex')
ylabel('$f_2(\chi)$','interpreter','latex')
title('Pareto front')
box on
% axis([15 50 -0.01 1.5])
legend('Nondominated points','Sensitivity analysis','Best for model')


%% sensitivity
mean(sensitiveF(:,1)./fval(sensitiveID,1))-1
mean(sensitiveF(:,2)./fval(sensitiveID,2))-1

%% goodness of fitting, rmse
sqrt(fval(sensitiveID,1))
sqrt(fval(sensitiveID,2).^2)





