clc; clear; close all;
 munlock('UWerr') 
 clear('UWerr') 
 [traces, variances, num_configs, num_timeslices, config_numbers] = read_raw;

% Compute correlator C(tau) and noise variance
max_tau = 4;
C_tau = zeros(max_tau+1, 1);
C_tau_err = zeros(max_tau+1, 1);

for tau = 0:max_tau
    val_vec = zeros(num_configs, 1);
    noise_var_vec = zeros(num_configs, 1);

    for idx = 1:num_configs
        tr = traces(:, idx);
        var = variances(:, idx)/16;

        % Periodic shift using circshift
        tr_shift = circshift(tr, -tau);
        var_shift = circshift(var, -tau);

        % Mean correlator for this configuration
        val = mean(tr .* tr_shift);

        % Stochastic variance (var_noise = 1/T^2 sum [...])
        var_noise = mean(var .* (tr_shift.^2) + var_shift .* (tr.^2)) / num_timeslices^2;
                
        var_deriv = 0;
        for tp = 1:num_timeslices
            tp_minus = mod(tp - tau - 1, num_timeslices) + 1;
            tp_plus = mod(tp + tau - 1, num_timeslices) + 1;
            coeff = (tr(tp_plus) + tr(tp_minus));% / num_timeslices;
            var_deriv = var_deriv + coeff^2 * var(tp);
        end

        var_deriv = var_deriv / num_timeslices^2;
        var_noise = var_deriv;
        val_vec(idx) = val;
        noise_var_vec(idx) = var_noise;
    end

    % Apply UWerr to get gauge error
    [mean_val, dval] = UWerr(val_vec, 1.5, [], [], 1);
     
    % Final total error: combine UWerr variance and noise estimate (Eq. 13)
    mean_noise_var = mean(noise_var_vec);
    total_var = dval^2 + mean_noise_var / (num_configs);

    C_tau(tau+1) = mean_val;
    C_tau_err(tau+1) = sqrt(total_var);
end

% Plot only C(tau) with error bars
figure;
errorbar(0:max_tau, C_tau, C_tau_err, 'o-', 'LineWidth', 1.2);
xlabel('$\\tau$', 'Interpreter', 'latex');
ylabel('$C(\\tau)$', 'Interpreter', 'latex');
title('Disconnected correlator $C(\\tau)$ with full variance', 'Interpreter', 'latex');
set(gca, 'YScale', 'log');
grid on;

% Save figure
savefig('correlator.fig');