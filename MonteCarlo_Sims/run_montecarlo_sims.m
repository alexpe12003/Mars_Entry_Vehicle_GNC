% =====================================================================
% Monte Carlo Simulation for Entry Vehicle Landing Dispersion Analysis
% =====================================================================
%
% This script performs a Monte Carlo simulation with N = 500 runs to assess
% the sensitivity of the entry vehicle's landing performance to variations
% in initial conditions and physical parameters. Each run perturbs key 
% variables such as initial velocity, flight path angle, mass, and attitude.
%
% The simulation workflow consists of:
%   1. Generating randomized initial conditions for each run
%   2. Assigning them to the base workspace for use in the Simulink model
%   3. Running the entry vehicle model (entryVehicle.slx)
%   4. Recording final landing coordinates and downrange performance
%   5. Plotting the dispersion of landing points in latitude (longitude fixed)
%   6. Visualizing the distribution of downrange errors with a histogram
%   7. Computing statistical metrics (mean and standard deviation of error)
%
% Outputs:
%   - Landing dispersion plot (latitude vs fixed longitude)
%   - Histogram of downrange errors (relative to nominal)
%   - Printed statistics on downrange performance
%
% Notes:
%   - The nominal target is defined by lat_c and lon_c
%   - Final landing states are extracted from the Simulink output
%   - Some parameters are commented out and can be re-enabled for further studies
%
% Author: [Your Name]
% Date: [Current Date]
% =====================================================================

N = 100;  % Number of Monte Carlo runs
results = struct();  % Store simulation results

% Nominal landing target
lat_c = -0.417;%-1.76;     % deg
lon_c = 15.495;%19.13;     % deg

lats = zeros(1, N);
lons = zeros(1, N);
s_star_nominal=1.75364*10^6; %2.209934992639319e+06;

% Initialize variables for best, worst success, and worst failure cases
best_case = inf;
worst_success_case = -inf;
worst_failure_case = inf;

best_case_s_final = 0;
worst_success_case_s_final = 0;
worst_failure_case_s_final = 0;

% Track the successful and failure cases for downrange error
successful_errors = [];
failure_errors = [];

for i = 1:N
    % Display the current iteration in the console
    fprintf('Running iteration %d of %d...\n', i, N);
    % === 1. Perturb initial conditions ===
    V_0_var       = 11e3 * (1 + 0.02*randn());
    gamma_0_var   = deg2rad(-9.536 + 0.1*randn());
    chi_0_var     = deg2rad(90 + 1.0*randn());  % heading
    R_0_var       = (220e3 + Mars_radius) + 3000*randn();
    tau_0_var     = deg2rad(0.1*randn());       % ±0.2 deg longitude error
    lambda_0_var  = deg2rad(0.1*randn());       % ±0.2 deg latitude error

    alpha_0_var   = deg2rad(-23.82 + 1.0*randn());
    beta_0_var    = deg2rad(0 + 0.5*randn());
    sigma_0_var   = deg2rad(100 + 5*randn());

    p_0_var       = deg2rad(0.2*randn());
    q_0_var       = deg2rad(0.2*randn());
    r_0_var       = deg2rad(0.2*randn());
    m_0_var       = 4976 * (1 + 0.05*randn());

    I_xx_var      = 5617.61 * (1 + 0.05*randn());
    I_yy_var      = 4454.62 * (1 + 0.05*randn());
    I_zz_var      = 4454.80 * (1 + 0.05*randn());

    % === 2. Assign to base workspace ===
    assignin('base', 'V_0', V_0_var);
    assignin('base', 'gamma_0', gamma_0_var);
    assignin('base', 'R_0', R_0_var);
    assignin('base', 'sigma_0', sigma_0_var);
    assignin('base', 'm_0', m_0_var);

    % === 3. Run Simulink model ===
    simOut = sim('entryVehicle.slx');
    state = simOut.state.signals.values;
    
    % === Extract final downrange value from simulation ===
    downrange = simOut.downrange.signals.values;  % should be [Nx1]
    s_final = downrange(end);  % final value in meters
    
    % Store downrange error relative to nominal 
    s_error_pct = 100 * (s_final - s_star_nominal) / s_star_nominal;  % error percentage
    s_errors_pct(i) = s_error_pct;
    results(i).s_final = s_final;
    results(i).s_error_pct = s_error_pct;

    % === Track the best, worst success, and worst failure cases ===
    if s_errors_pct(i) < 100
        % Successful landing (error < 100%)
        successful_errors = [successful_errors, s_errors_pct(i)];
        
        % Track the worst success case (largest error in successful landings)
        if s_errors_pct(i) > worst_success_case
            worst_success_case = s_errors_pct(i);
            worst_success_case_s_final = s_final;
        end
    else
        % Failure (error > 100%)
        failure_errors = [failure_errors, s_errors_pct(i)];
        
        % Track the worst failure case (largest error in failures)
        if s_errors_pct(i) < worst_failure_case
            worst_failure_case = s_errors_pct(i);
            worst_failure_case_s_final = s_final;
        end
    end

    % Track the best case (minimum error, can be a success or failure)
    if s_errors_pct(i) < best_case
        best_case = s_errors_pct(i);
        best_case_s_final = s_final;
    end

    % === 4. Extract final longitude and latitude ===
    final_tau = state(end, 2);      % Longitude [rad]
    final_delta = state(end, 3);    % Latitude [rad]

    lons(i) = rad2deg(final_tau);
    lats(i) = rad2deg(final_delta);

    % Optional results
    results(i).final_lat = lats(i);
    results(i).final_lon = lons(i);
end

% === 5. Plot Landing Dispersion (Longitude Fixed) ===
figure;
hold on; grid on;

% Force all longitudes to the nominal value for plotting
lons = lon_c * ones(size(lons));  % Override longitudes

% Plot Monte Carlo landing points
plot(lons, lats, 'b.', 'MarkerSize', 10, 'DisplayName', 'Landing Points');

% Plot centroid (mean landing)
lat_mean = mean(lats);
lon_mean = lon_c;
plot(lon_mean, lat_mean, 'ko', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Mean Landing');

% Plot nominal target point
plot(lon_c, lat_c, 'rx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Nominal Target');

% Draw 50 km radius circle around nominal point
r_nominal = 50000;  % 50 km in meters
deg_per_m_lat = 1 / 111320;
deg_per_m_lon = 1 / (111320 * cosd(lat_c));

theta = linspace(0, 2*pi, 300);
lat_circle_nom = lat_c + r_nominal * deg_per_m_lat * sin(theta);
lon_circle_nom = lon_c + r_nominal * deg_per_m_lon * cos(theta);
plot(lon_circle_nom, lat_circle_nom, 'r--', 'LineWidth', 2, 'DisplayName', '100 km Nominal Circle');

% Final plot styling
title('Landing Dispersion with Bank Angle Modulation');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Location', 'bestoutside');
axis equal;

% === 6. Plot Histogram ===
figure;
histogram(abs(s_errors_pct), 'BinWidth', 1, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Downrange Error [%]');
ylabel('Frequency');
title('Histogram of Downrange Error (% of Nominal)');
grid on;

% === 7. Print Valuable Results ===
fprintf('Best Case: %.2f%% error | Final Landing Distance = %.2f m\n', best_case, best_case_s_final);
fprintf('Worst Success Case: %.2f%% error | Final Landing Distance = %.2f m\n', worst_success_case, worst_success_case_s_final);
fprintf('Worst Failure Case: %.2f%% error | Final Landing Distance = %.2f m\n', worst_failure_case, worst_failure_case_s_final);

% Successful landing statistics
if ~isempty(successful_errors)
    mean_success_err = mean(successful_errors);
    std_success_err = std(successful_errors);
    fprintf('Successful Landing Stats - Mean = %.2f%% | Std = %.2f%%\n', mean_success_err, std_success_err);
else
    fprintf('No successful landings (errors < 100)\n');
end

% Failure landing statistics
if ~isempty(failure_errors)
    mean_failure_err = mean(failure_errors);
    std_failure_err = std(failure_errors);
    fprintf('Failure Landing Stats - Mean = %.2f%% | Std = %.2f%%\n', mean_failure_err, std_failure_err);
else
    fprintf('No failures (errors > 100)\n');
end
