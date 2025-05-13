close all

% Compute test and test_var
tr = traces(:,1);
tr_err = sqrt(variances(:,1) / 16); 

% X-axis values
x = 1:length(tr);  

% Plot with red error bars and filled black circles
h = errorbar(x, tr, tr_err, 'o', ...
    'DisplayName', '$\mathrm{tr}((\gamma_5 D(t,t))^{-1})$', ...
    'LineWidth', 1.2, ...
    'Color', 'r', ...           
    'MarkerEdgeColor', 'k', ... 
    'MarkerFaceColor', 'k');    

legend show
legend('Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex')
ylabel('$\mathrm{tr}((\gamma_5 D(t,t))^{-1})$', 'Interpreter', 'latex')
title('Traces at each time slice $t$', 'Interpreter', 'latex')
grid on
