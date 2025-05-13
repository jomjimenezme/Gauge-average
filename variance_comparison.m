values = zeros(num_configs, 1);
errors_corr = zeros(num_configs, 1);
errors_derivative = zeros(num_configs, 1);  % Derivative-based method

for idx = 1:num_configs
    tr = traces(:, idx);
    var = variances(:, idx)/16;

    tr_shift = circshift(tr, -tau);
    var_shift = circshift(var, -tau);

    mu_t = tr;
    mu_tp = tr_shift;
    sigma2_t = var;
    sigma2_tp = var_shift;

    % Estimate covariance directly
    cov_tp = zeros(num_timeslices, 1);
    for t = 1:num_timeslices
        t_shift = mod(t + tau - 1, num_timeslices) + 1;
        cov_tp(t) = (tr(t) * tr(t_shift)) - mu_t(t) * mu_tp(t);
    end

    % Apply correlated variance formula
    var_corr = mean(mu_t.^2 .* sigma2_tp + ...
                    mu_tp.^2 .* sigma2_t + ...
                    sigma2_t .* sigma2_tp + ...
                    2 * mu_t .* mu_tp .* cov_tp) / num_timeslices^2;

    % Derivative-based method 
    var_deriv = 0;
    for tp = 1:num_timeslices
        tp_minus = mod(tp - tau - 1, num_timeslices) + 1;
        tp_plus = mod(tp + tau - 1, num_timeslices) + 1;
        coeff = (tr(tp_plus) + tr(tp_minus));
        var_deriv = var_deriv + coeff^2 * var(tp);
    end
    var_deriv = var_deriv / num_timeslices^2;

    values(idx) = mean(tr .* tr_shift);
    errors_corr(idx) = sqrt(var_corr);
    errors_derivative(idx) = sqrt(var_deriv);
end

figure;
hold on;
errorbar(config_numbers, values, errors_corr, 's--', 'LineWidth', 1.2, 'DisplayName', 'Correlated product variance');
errorbar(config_numbers, values, errors_derivative, 'd-.', 'LineWidth', 1.2, 'DisplayName', 'Derivative formula');
xlabel('Configuration number', 'Interpreter', 'latex');
ylabel('Correlator at $$\\tau = 0$$', 'Interpreter', 'latex');
title('Comparison of stochastic error estimates at fixed $$\\tau$$', 'Interpreter', 'latex');
legend('Location', 'best');
grid on;

savefig('val_vs_config_comparison.fig');
