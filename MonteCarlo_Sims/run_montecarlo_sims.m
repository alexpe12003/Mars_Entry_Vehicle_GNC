N = 10;  % Number of Monte Carlo runs
results = struct();  % Store simulation results

% Nominal landing target
lat_c = -1.76;     % deg
lon_c = 19.13;     % deg

lats = zeros(1, N);
lons = zeros(1, N);
s_star_nominal=2.209934992639319e+06;
for i = 1:N
    % === 1. Perturb initial conditions ===
    V_0_var       = 11e3 * (1 + 0.02*randn());
    gamma_0_var   = deg2rad(-9.536 + 0.1*randn());
    chi_0_var     = deg2rad(90 + 1.0*randn());  % heading
    R_0_var       = (220e3 + Mars_radius) + 3000*randn();
    %tau_0_var     = deg2rad(0.2*randn());       % ±0.2 deg longitude error
    %lambda_0_var  = deg2rad(0.2*randn());       % ±0.2 deg latitude error

    %alpha_0_var   = deg2rad(-23.82 + 1.0*randn());
    %beta_0_var    = deg2rad(0 + 0.5*randn());
    sigma_0_var   = deg2rad(100 + 5*randn());

    %p_0_var       = deg2rad(0.2*randn());
    %q_0_var       = deg2rad(0.2*randn());
    %r_0_var       = deg2rad(0.2*randn());

    %I_xx_var      = 5617.61 * (1 + 0.05*randn());
    %I_yy_var      = 4454.62 * (1 + 0.05*randn());
    %I_zz_var      = 4454.80 * (1 + 0.05*randn());

    %I_xy_var      = 0;
    %I_yz_var      = 0;
    %I_xz_var      = 0;

    %S_ref_var     = 12 * (1 + 0.03*randn());
    %d_ref_var     = 3.9 * (1 + 0.03*randn());

    %r_cm_var      = [-0.137; 0; 1.8] + 0.01*randn(3,1);  % up to 1 cm deviation

    % === 2. Assign to base workspace ===
    assignin('base', 'V_0', V_0_var);
    assignin('base', 'gamma_0', gamma_0_var);
    %assignin('base', 'chi_0', chi_0_var);
    assignin('base', 'R_0', R_0_var);
    %assignin('base', 'tau_0', tau_0_var);
    %assignin('base', 'lambda_0', lambda_0_var);

    %assignin('base', 'alpha_0', alpha_0_var);
    %assignin('base', 'beta_0', beta_0_var);
    assignin('base', 'sigma_0', sigma_0_var);

    %assignin('base', 'p_0', p_0_var);
    %assignin('base', 'q_0', q_0_var);
    %assignin('base', 'r_0', r_0_var);

    %assignin('base', 'm_0', m_0_var);
    %assignin('base', 'I_xx', I_xx_var);
    %assignin('base', 'I_yy', I_yy_var);
    %assignin('base', 'I_zz', I_zz_var);
    %assignin('base', 'I_xy', I_xy_var);
    %assignin('base', 'I_yz', I_yz_var);
    %assignin('base', 'I_xz', I_xz_var);

    %assignin('base', 'S_ref', S_ref_var);
    %assignin('base', 'd_ref', d_ref_var);
    %assignin('base', 'r_cm', r_cm_var);

    % === 3. Run Simulink model ===
    simOut = sim('entryVehicle.slx');
    state = simOut.state.signals.values;
    % === Extract final downrange value from simulation ===
    downrange = simOut.downrange.signals.values;  % should be [Nx1]
    s_final = downrange(end);  % final value in meters
    
    % Store downrange error relative to nominal
    s_error_pct = 100 * (s_final - s_star_nominal) / s_star_nominal;  % erro percentual
    s_errors_pct(i) = s_error_pct;
    results(i).s_final = s_final;
    results(i).s_error_pct = s_error_pct;

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

% === Draw 100 km radius circle around nominal point ===
r_nominal = 100000;  % 10 km in meters
deg_per_m_lat = 1 / 111320;
deg_per_m_lon = 1 / (111320 * cosd(lat_c));

theta = linspace(0, 2*pi, 300);
lat_circle_nom = lat_c + r_nominal * deg_per_m_lat * sin(theta);
lon_circle_nom = lon_c + r_nominal * deg_per_m_lon * cos(theta);
plot(lon_circle_nom, lat_circle_nom, 'r--', 'LineWidth', 2, 'DisplayName', '100 km Nominal Circle');

% === Final plot styling ===
title('Landing Dispersion with Bank Angle Modulation');
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Location', 'bestoutside');
axis equal;

% === 4. Plotar histograma ===
figure;
histogram(abs(s_errors_pct), 'BinWidth', 1, 'FaceColor', [0.2 0.4 0.8]);
xlabel('Downrange Error [%]');
ylabel('Frequency');
title('Histogram of Downrange Error (% of Nominal)');
grid on;

% === 5. Estatísticas ===
mean_err = mean(abs(s_errors_pct));
std_err = std(abs(s_errors_pct));
fprintf('Downrange Error: Mean = %.2f%% | Std = %.2f%%\n', mean_err, std_err);