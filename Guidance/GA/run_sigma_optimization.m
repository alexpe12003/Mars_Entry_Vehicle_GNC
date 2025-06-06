function run_sigma_optimization()
    N = 10;                            % Number of control points (sigma values)
    sigma_bounds = deg2rad(100);       % [-60°, 60°]
    lb = -sigma_bounds * ones(1, N);
    ub =  sigma_bounds * ones(1, N);

    % Genetic Algorithm
    options = optimoptions('ga', ...
        'PopulationSize', 30, ...
        'MaxGenerations', 50, ...
        'UseParallel', false, ...
        'Display', 'iter');

    [best_sigma, best_heat] = ga(@(x) simulate_qc_cost(x), N, [], [], [], [], lb, ub, [], options);

    fprintf("Best total heat flux: %.2f W/m²\n", best_heat);
    disp("Best sigma reference (rad):")
    disp(best_sigma)

    % Save best reference as a timeseries or signal
    t_total = 500;
    t_ref = linspace(0, t_total, N);
    sigma_signal = timeseries(best_sigma, t_ref);
    save('sigmaDataset.mat', 'sigma_signal');
end


