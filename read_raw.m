function [traces, variances, num_configs, num_timeslices, config_numbers] = read_raw

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



end