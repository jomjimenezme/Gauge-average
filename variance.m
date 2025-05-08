clc; clear; close all;

% List of configuration numbers
config_numbers = [20, 40, 60, 80, 120, 140, 160, 180, 200, ...
                  220, 240, 260, 280, 300, 320, 340, 360, 380, ...
                  420, 440, 460, 480];

num_configs = numel(config_numbers);
num_timeslices = 192;
num_noise = 16;

% Initialize traces and variances
traces = zeros(num_timeslices, num_configs);
variances = zeros(num_timeslices, num_configs);

% Read all traces and variances
for idx = 1:num_configs
    config = config_numbers(idx);
    file_name = sprintf('%d/MLMC_latest/mpi-outsh_%d', config, config);

    filetext = fileread(file_name);
    fprintf("%s\n", file_name);

    expr_trace = 'trace: ([-+]?\d+(\.\d+)?)';
    matches_trace = regexp(filetext, expr_trace, 'tokens');
    matches_trace = cellfun(@(x) x{1}, matches_trace, 'UniformOutput', false);
    traces(:, idx) = str2double(matches_trace);

    expr_var = 'variance: ([-+]?\d+(\.\d+)?)';
    matches_var = regexp(filetext, expr_var, 'tokens');
    matches_var = cellfun(@(x) x{1}, matches_var, 'UniformOutput', false);
    variances(:, idx) = str2double(matches_var);
end

% Compute correlator C(tau) and noise variance
max_tau = 7;
C_tau = zeros(max_tau+1, 1);
C_tau_err = zeros(max_tau+1, 1);

for tau = 0:max_tau
    val_vec = zeros(num_configs, 1);
    noise_var_vec = zeros(num_configs, 1);

    for idx = 1:num_configs
        tr = traces(:, idx);
        var = variances(:, idx);

        % Periodic shift using circshift
        tr_shift = circshift(tr, -tau);
        var_shift = circshift(var, -tau);

        % Mean correlator for this configuration
        val = mean(tr .* tr_shift);

        % Stochastic variance (Eq. var_noise = 1/T^2 sum [...])
        var_noise = mean(var .* (tr_shift.^2) + var_shift .* (tr.^2)) / num_timeslices^2;

        val_vec(idx) = val;
        noise_var_vec(idx) = var_noise;
    end

    % Apply UWerr to get gauge error
    [mean_val, dval] = UWerr(val_vec, 1.5, [], [], 1);

    % Final total error: combine UWerr variance and noise estimate (Eq. 13)
    mean_noise_var = mean(noise_var_vec);
    total_var = dval^2 + mean_noise_var / (num_configs * num_noise);

    C_tau(tau+1) = mean_val;
    C_tau_err(tau+1) = sqrt(total_var);
end

% Plot only C(tau) with error bars
figure;
errorbar(0:max_tau, C_tau, C_tau_err, 'o-', 'LineWidth', 1.2);
xlabel('$\\tau$', 'Interpreter', 'latex');
ylabel('$C(\\tau)$', 'Interpreter', 'latex');
title('Disconnected correlator $C(\\tau)$ with full variance', 'Interpreter', 'latex');
grid on;

% Save figure
savefig('correlator.fig');