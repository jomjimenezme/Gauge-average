clc; clear; close all;

% List of configuration numbers
config_numbers = [20, 40, 60, 80, 120, 140, 160, 180, 200, ...
                  220, 240, 260, 280, 300, 320, 340, 360, 380, ...
                  420, 440, 460, 480];

num_configs = numel(config_numbers);
num_timeslices = 192;

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

% Choose fixed tau (e.g. tau = 0)
tau = 0;

values = zeros(num_configs, 1);
errors_taylor = zeros(num_configs, 1);
errors_uncorr = zeros(num_configs, 1);

for idx = 1:num_configs
    tr = traces(:, idx);
    var = variances(:, idx);

    tr_shift = circshift(tr, -tau);
    var_shift = circshift(var, -tau);

    val = mean(tr .* tr_shift);
    
    % Taylor expansion (2-term)
    var_taylor = mean(var .* (tr_shift.^2) + var_shift .* (tr.^2)) / num_timeslices^2;

    % Full uncorrelated product variance (3-term)
    var_uncorr = mean(var .* (tr_shift.^2) + var_shift .* (tr.^2) + var .* var_shift) / num_timeslices^2;

    values(idx) = val;
    errors_taylor(idx) = sqrt(var_taylor);
    errors_uncorr(idx) = sqrt(var_uncorr);
end

% Plot value vs config number with both error bars
figure;
hold on;
errorbar(config_numbers, values, errors_taylor, 'o-', 'LineWidth', 1.2, 'DisplayName', 'Taylor ');
errorbar(config_numbers, values, errors_uncorr, 's--', 'LineWidth', 1.2, 'DisplayName', 'Uncorrelated product ');
xlabel('Configuration number', 'Interpreter', 'latex');
ylabel('Correlator at $$\tau = 0$$', 'Interpreter', 'latex');
title('Comparison of stochastic error estimates at fixed $$\tau$$', 'Interpreter', 'latex');
legend('Location', 'best');
grid on;

% Save figure
savefig('val_vs_config_comparison.fig');
