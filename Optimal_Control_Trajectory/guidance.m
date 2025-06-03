simOut = sim("entryVehicle.slx");

m_0 = 4976;
Sref = 12;
Mars_radius = 6371e3;

% === Constants for new gain computation ===
ballistic_beta = 340;  % [kg/m^2]
L_over_D = 0.3;

% === Extract state vectors ===
state = simOut.state.signals.values;
t     = simOut.state.time;

v_star     = state(:, 7);
gamma_star = state(:, 8);
sigma      = state(:, 6);
h_star     = state(:, 1) - Mars_radius;
hdot_star  = v_star .* sin(gamma_star);
s_dot_star = v_star .* cos(gamma_star);

% === Integrate s_dot using ODE solver for higher accuracy ===
s_dot_fun = @(tt) interp1(t, s_dot_star, tt, 'linear', 'extrap');
[~, s_sol] = ode45(@(tt, ss) s_dot_fun(tt), t, 0);
s_star = s_sol;

% === Extract drag and atmospheric data ===
DLS     = simOut.D_L_S.signals.values;
D_star  = DLS(:, 1);
rho     = simOut.rho.signals.values;
gs      = simOut.gs.signals.values / m_0;

% === Compute gains from costate equations (updated version) ===
gain_table = compute_F1_F2_F3_from_costates( ...
    v_star, gamma_star, h_star, hdot_star, D_star, ...
    Mars_radius, m_0, rho, gs, Sref, sigma, ballistic_beta, L_over_D);

% === Filter for velocity range of interest ===
v_min = 0;
v_max = 11000;
valid_idx = (v_star >= v_min) & (v_star <= v_max);

% Apply filtering
v_star     = v_star(valid_idx);
hdot_star  = hdot_star(valid_idx);
s_star     = s_star(valid_idx);
D_star     = D_star(valid_idx);
F1         = gain_table.F1(valid_idx);
F2         = gain_table.F2(valid_idx);
F3         = gain_table.F3(valid_idx);
t          = t(valid_idx);

% === Bin and average by velocity ===
v_bin_width = 10;
v_binned = round(v_star / v_bin_width) * v_bin_width;
[v_unique_ds, ~, bin_idx] = unique(v_binned);

hdot_lut_ds = accumarray(bin_idx, hdot_star, [], @mean);
s_lut_ds    = accumarray(bin_idx, s_star,    [], @mean);
drag_lut_ds = accumarray(bin_idx, D_star,    [], @mean);
F1_ds       = accumarray(bin_idx, F1, [], @mean);
F2_ds       = accumarray(bin_idx, F2, [], @mean);
F3_ds       = accumarray(bin_idx, F3, [], @mean);

% === Pack results into gain table structure ===
gain_table_ds.v  = v_unique_ds;
gain_table_ds.F1 = F1_ds;
gain_table_ds.F2 = F2_ds;
gain_table_ds.F3 = F3_ds;

save('guidance_lookup_tables.mat', ...
    'v_unique_ds', 'hdot_lut_ds', 's_lut_ds', 'drag_lut_ds', 'gain_table_ds');

% === Sanity check plots ===
figure('Name','Guidance Lookup Tables (Raw, Binned)', 'NumberTitle','off');

subplot(2,3,1)
plot(v_unique_ds, hdot_lut_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('hdot^* [m/s]');
title('Altitude Rate vs Velocity');

subplot(2,3,2)
plot(v_unique_ds, s_lut_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('s^* [m]');
title('Downrange vs Velocity');

subplot(2,3,3)
plot(v_unique_ds, drag_lut_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('D^* [N]');
title('Drag vs Velocity');

subplot(2,3,4)
plot(v_unique_ds, F1_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F1');
title('Gain F1 vs Velocity');

subplot(2,3,5)
plot(v_unique_ds, F2_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F2');
title('Gain F2 vs Velocity');

subplot(2,3,6)
plot(v_unique_ds, F3_ds, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F3');
title('Gain F3 vs Velocity');

% === Downrange vs Time ===
figure('Name','Downrange vs Time','NumberTitle','off');
plot(t, s_star, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('s^* [m]');
title('Downrange vs Time');

