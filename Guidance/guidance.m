% =====================================================================
% Guidance Gain Computation and Lookup Table Generation
% =====================================================================
%
% This script runs the open-loop entry vehicle simulation and processes
% the resulting trajectory data to compute the costate-based guidance gains
% (F1, F2, F3). These gains are binned by velocity to build a set of
% lookup tables that can later be used for real-time entry guidance.
%
% Steps performed:
%   1. Load planetary and atmospheric constants for Earth.
%   2. Run the entry vehicle Simulink model and extract state data.
%   3. Compute vertical and horizontal velocity components (hdot, sdot).
%   4. Integrate sdot to obtain downrange position s(t).
%   5. Extract drag and density profiles from the simulation.
%   6. Call a function to compute gains F1, F2, F3 from costate equations.
%   7. Filter data within a valid velocity range and bin it.
%   8. Average each variable within velocity bins to smooth results.
%   9. Save the guidance lookup tables for later use in the controller.
%  10. Generate sanity check plots to validate table behavior.
%
% Notes:
%   - The simulation uses Earth parameters but reuses Mars-style notation
%     (e.g., `Mars_radius`) for compatibility with prior models.
%   - Final results are saved to the 'Guidance/data_sets_guidance/' folder.
%
% Authors:
%   Alexandre Alves Pereira       Nº112288
%   Eduardo De Almeida Helena     Nº102793
%   Guilherme Alexandre Pires Martins Nº112118
%   Rafael Coimbra Azeiteiro      Nº102478
%   Rodrigo Marques Pereira       Nº112117
%
% Professor:
%   Alain de Souza
%
% Course Project — Guidance, Navigation and Control, IST Lisbon University
% =====================================================================
%%
%To run another trajectory please uncomment the code.
%% Planet Constants - Earth
% 
% % === Earth Radius ===
 Mars_radius = 6371e3;  % Radius of Earth in meters [m]
 H_s = 7050;  % Scale height of Earth' atmosphere [m]
 rho_0 = 1.225;  % Reference density of Earth' atmosphere at sea level [kg/m^3]
% 
% %% Simulation Setup
% 
% simOut = sim("entryVehicle.slx");
% %sample_rate_guidance
% % === Constants for New Gain Computation ===
m_0 = 4976;  % [kg] - Initial mass
Sref = 12;   % [m^2] - Reference area
% 
ballistic_beta = 340;  % [kg/m^2] - Ballistic coefficient
L_over_D = 0.3;       % [Dimensionless] - Lift-to-Drag ratio
dynamics_constants = [m_0, S_ref, H_s, ballistic_beta, rho_0, L_over_D, Mars_radius];
% 
% % === Extract State Vectors from Simulation ===
% state = simOut.state.signals.values;
% t = simOut.state.time;
% 
% v_star = state(:, 7);  % Velocity
% gamma_star = state(:, 8);  % Flight-path angle
% sigma = state(:, 6);  % Bank angle
% h_star = state(:, 1) - Mars_radius;  % Altitude relative to Mars' surface
% hdot_star = v_star .* sin(gamma_star);  % Vertical velocity component
% s_dot_star = v_star .* cos(gamma_star);  % Horizontal velocity component
% 
% % === Integrate s_dot Using ODE Solver for Higher Accuracy ===
% s_dot_fun = @(tt) interp1(t, s_dot_star, tt, 'linear', 'extrap');
% [~, s_sol] = ode45(@(tt, ss) s_dot_fun(tt), t, 0);
% s_star = s_sol;  % Downrange position
% 
% % === Extract Drag and Atmospheric Data from Simulation ===
% DLS = simOut.D_L_S.signals.values;
% D_star = DLS(:, 1);  % Drag force [N]
% rho = simOut.rho.signals.values;  % Atmospheric density [kg/m^3]
% gs = simOut.gs.signals.values / m_0;  % Gravitational acceleration per unit mass
% 
% % === Compute Gains from Costate Equations ===
% gain_table = compute_F1_F2_F3_from_costates( ...
%     v_star, gamma_star, h_star, hdot_star, D_star, ...
%     Mars_radius, m_0, rho, gs, Sref, sigma, ballistic_beta, L_over_D,t);
% 
% % === Filter for Velocity Range of Interest ===
% v_min = 1000;  % Minimum velocity [m/s]
% v_max = 10000;  % Maximum velocity [m/s]
% valid_idx = (v_star >= v_min) & (v_star <= v_max);  % Filter condition
% 
% % Apply filtering to all relevant variables
% v_star = v_star(valid_idx);
% hdot_star = hdot_star(valid_idx);
% s_star = s_star(valid_idx);
% D_star = D_star(valid_idx);
% F1 = gain_table.F1(valid_idx);
% F2 = gain_table.F2(valid_idx);
% F3 = gain_table.F3(valid_idx);
% t = t(valid_idx);
% sigma=sigma(valid_idx);
% h_star=h_star(valid_idx);
% 
% % === Bin and Average by Velocity ===
% v_bin_width = 10;  % Velocity bin width [m/s]
% v_binned = round(v_star / v_bin_width) * v_bin_width;  % Binning velocity values
% [v_unique_ds, ~, bin_idx] = unique(v_binned);  % Unique velocity bins
% 
% % Bin data and compute averages
% hdot_lut_ds = accumarray(bin_idx, hdot_star, [], @mean);
% s_lut_ds = accumarray(bin_idx, s_star, [], @mean);
% drag_lut_ds = accumarray(bin_idx, D_star, [], @mean);
% F1_ds = accumarray(bin_idx, F1, [], @mean);
% F2_ds = accumarray(bin_idx, F2, [], @mean);
% F3_ds = accumarray(bin_idx, F3, [], @mean);
% 
% % === Pack Results into Gain Table Structure ===
% gain_table_ds.v = v_unique_ds;
% gain_table_ds.F1 = F1_ds;
% gain_table_ds.F2 = F2_ds;
% gain_table_ds.F3 = F3_ds;
% 
% % Now save the data
% save('Guidance/data_sets_guidance/guidance_lookup_tables.mat', 'v_unique_ds', 'hdot_lut_ds', 'drag_lut_ds','s_lut_ds', 'gain_table_ds');
% 
% %% Sanity Check Plots
% 
% % Plot raw and binned guidance lookup tables
% figure('Name','Guidance Lookup Tables (Raw, Binned)', 'NumberTitle','off');
% 
% subplot(2,3,1)
% plot(v_unique_ds, hdot_lut_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('hdot^* [m/s]');
% title('Altitude Rate vs Velocity');
% 
% subplot(2,3,2)
% plot(v_unique_ds, s_lut_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('s^* [m]');
% title('Downrange vs Velocity');
% 
% subplot(2,3,3)
% plot(v_unique_ds, drag_lut_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('D^* [N]');
% title('Drag vs Velocity');
% 
% subplot(2,3,4)
% plot(v_unique_ds, F1_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('F1');
% title('Gain F1 vs Velocity');
% 
% subplot(2,3,5)
% plot(v_unique_ds, F2_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('F2');
% title('Gain F2 vs Velocity');
% 
% subplot(2,3,6)
% plot(v_unique_ds, F3_ds, 'LineWidth', 1.2); grid on;
% xlabel('Velocity [m/s]'); ylabel('F3');
% title('Gain F3 vs Velocity');
% 
% % === Downrange vs Time ===
% figure('Name','Downrange vs Time','NumberTitle','off');
% plot(t, s_star, 'LineWidth', 1.5);
% xlabel('Time [s]');
% ylabel('s^* [m]');
% title('Downrange vs Time');
% 
% % === Bank Angle (sigma) vs Time ===
% figure('Name','Bank Angle (sigma) vs Time','NumberTitle','off');
% plot(t, rad2deg(sigma), 'LineWidth', 1.5);
% xlabel('Time [s]');
% ylabel('Bank Angle (sigma) [deg]');
% title('Bank Angle (sigma) vs Time');
% grid on;
% 
% % === Altitude vs Downrange ===
% figure('Name','Downrange vs Time','NumberTitle','off');
% plot(s_star*10^-3, h_star*10^-3, 'LineWidth', 1.5);
% xlabel('Downrange[km]');
% ylabel('Altiude [km]');
% title('Altitude vs Downrange');
% 
% % === Bank Angle (sigma) vs Time ===
% figure('Name','Drag Profile vs Time','NumberTitle','off');
% plot(t, D_star, 'LineWidth', 1.5);
% xlabel('Time [s]');
% ylabel('Drag Profile [N]');
% title('Drag vs Time');
% grid on;

disp("✔️ Guidance lookup tables with pre-defined nominal trajectory loaded.");

