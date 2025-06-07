%% Planet Constants - Earth

% === Earth Radius ===
Mars_radius = 6371e3;  % Radius of Earth in meters [m]

%% Simulation Setup

simOut = sim("entryVehicle.slx");
sample_rate_guidance = 0.01;

% === Constants for New Gain Computation ===
m_0 = 4976;  % [kg] - Initial mass
Sref = 12;   % [m^2] - Reference area

ballistic_beta = 340;  % [kg/m^2] - Ballistic coefficient
L_over_D = 0.3;       % [Dimensionless] - Lift-to-Drag ratio
dynamics_constants = [m_0, Sref, H_s, ballistic_beta, rho_0, L_over_D, Mars_radius];

% === Extract State Vectors from Simulation ===
state = simOut.state.signals.values;
t = simOut.state.time;

v_star = state(:, 7);  % Velocity
gamma_star = state(:, 8);  % Flight-path angle
sigma = state(:, 6);  % Bank angle
h_star = state(:, 1) - Mars_radius;  % Altitude relative to Mars' surface
hdot_star = v_star .* sin(gamma_star);  % Vertical velocity component
s_dot_star = v_star .* cos(gamma_star);  % Horizontal velocity component

% === Integrate s_dot Using ODE Solver for Higher Accuracy ===
s_dot_fun = @(tt) interp1(t, s_dot_star, tt, 'linear', 'extrap');
[~, s_sol] = ode45(@(tt, ss) s_dot_fun(tt), t, 0);
s_star = s_sol;  % Downrange position

% === Extract Drag and Atmospheric Data from Simulation ===
DLS = simOut.D_L_S.signals.values;
D_star = DLS(:, 1);  % Drag force [N]
rho = simOut.rho.signals.values;  % Atmospheric density [kg/m^3]
gs = simOut.gs.signals.values / m_0;  % Gravitational acceleration per unit mass

% === Compute Gains from Costate Equations ===
gain_table = compute_F1_F2_F3_from_costates( ...
    v_star, gamma_star, h_star, hdot_star, D_star, ...
    Mars_radius, m_0, rho, gs, Sref, sigma, ballistic_beta, L_over_D, t);


% === Filter for Velocity Range of Interest ===
v_min = 300;  % Minimum velocity [m/s]
v_max = 8000;  % Maximum velocity [m/s]
valid_idx = (v_star >= v_min) & (v_star <= v_max);  % Filter condition

% Apply filtering to all relevant variables
v_star = v_star(valid_idx);
hdot_star = hdot_star(valid_idx);
s_star = s_star(valid_idx);
D_star = D_star(valid_idx);
F1 = gain_table.F1(valid_idx);
F2 = gain_table.F2(valid_idx);
F3 = gain_table.F3(valid_idx);
t = t(valid_idx);
sigma=sigma(valid_idx);

% === Pack Results into Gain Table Structure ===
gain_table_ds.v = v_star;
gain_table_ds.F1 = F1;
gain_table_ds.F2 = F2;
gain_table_ds.F3 = F3;

% Now save the data without any filtering or binning
save('Guidance/data_sets_guidance/guidance_lookup_tables.mat', 'v_star', 'hdot_star', 'D_star', 'gain_table_ds', 's_star');

%% Sanity Check Plots

% Plot raw and binned guidance lookup tables
figure('Name','Guidance Lookup Tables (Raw, Binned)', 'NumberTitle','off');

subplot(2,3,1)
plot(v_star, hdot_star, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('hdot^* [m/s]');
title('Altitude Rate vs Velocity');

subplot(2,3,2)
plot(v_star, s_star, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('s^* [m]');
title('Downrange vs Velocity');

subplot(2,3,3)
plot(v_star, D_star, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('D^* [N]');
title('Drag vs Velocity');

subplot(2,3,4)
plot(v_star, F1, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F1');
title('Gain F1 vs Velocity');

subplot(2,3,5)
plot(v_star, F2, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F2');
title('Gain F2 vs Velocity');

subplot(2,3,6)
plot(v_star, F3, 'LineWidth', 1.2); grid on;
xlabel('Velocity [m/s]'); ylabel('F3');
title('Gain F3 vs Velocity');

% === Downrange vs Time ===
figure('Name','Downrange vs Time','NumberTitle','off');
plot(t, s_star, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('s^* [m]');
title('Downrange vs Time');

% === Bank Angle (sigma) vs Time ===
figure('Name','Bank Angle (sigma) vs Time','NumberTitle','off');
plot(t, sigma, 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Bank Angle (sigma) [rad]');
title('Bank Angle (sigma) vs Time');
grid on;
