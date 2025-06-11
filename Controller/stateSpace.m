% =====================================================================
% Attitude Control Design and Gain Scheduling for Entry Vehicle
% =====================================================================
%
% This script handles the attitude controller design and gain scheduling
% for an entry vehicle using linearized models at selected trajectory times.
% It performs the following steps:
%
% 1. Loads cleaned open-loop trajectory simulation results.
% 2. Extracts key state variables and aerodynamic parameters.
% 3. (Optional) Linearizes the system at a selected time index (e.g., 100s).
% 4. Designs an LQR controller for the attitude subsystem.
% 5. Builds decoupled longitudinal and lateral controllers with constraints.
% 6. Performs gain scheduling by computing controller gains at multiple 
%    points along the trajectory.
% 7. Saves gain lookup tables (LUTs) for use in Simulink implementation.
% 8. Creates a sigma reference signal for simulation via the Signal Editor.
%
% Notes:
%   - The gain schedule is stored in 'Controller/data_sets_control/'.
%   - Linearization is based on the capsule properties and trajectory data.
%   - Many controller design steps are currently commented for modular testing.
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

clc;

% === Load Data ===
load("Trajectory_openLoop_cleaned.mat");

% === Extract Relevant Signals from the Output ===
state = out.state.signals.values;
D_L_S = out.D_L_S.signals.values;
F_g = out.F_g.signals.values;
mach = out.mach.signals.values;
p_dyna = out.q.signals.values;

% === Linearization Test for Trajectory at 100s ===
time_step = 100;  % [s]
idx = time_step / stepTime;

% === Define Capsule Properties ===
capsule.m0 = m_0;
capsule.Ixx = I_xx;
capsule.Iyy = I_yy;
capsule.Izz = I_zz;
capsule.Ixy = I_xy;
capsule.Iyz = I_yz;
capsule.Ixz = I_xz;
capsule.S_ref = S_ref;
capsule.d_ref = d_ref;

% === Linearization of System for 100s ===
%[A, B, ~, ~] = capsule_linearize(capsule, idx, state, D_L_S, F_g, mach, p_dyna, rho_0, H_s, Mars_radius);

% === Extract Attitude Dynamics Subsystem ===
%A_attitude = A(4:9, 4:9);
%B_attitude = B(4:9, :);

% === LQR Controller Design ===
% State constraints
%pmax = deg2rad(3);    
%qmax = deg2rad(3);    
%rmax = deg2rad(2);    
%alfamax = deg2rad(3); 
%betamax = deg2rad(3); 
%sigmamax = deg2rad(3);

% Actuator limits (moments)
%mx_max = 1000; 
%my_max = 1000; 
%mz_max = 1000; 

% Build Q and R matrices for LQR design
%Q = diag([ ...
 %   1/(pmax)^2, ...
 %   1/(qmax)^2, ...
 %   1/(rmax)^2, ...
 %   1/(alfamax)^2, ...
 %   1/(betamax)^2, ...
 %   1/(sigmamax)^2]);

%R = diag([ ...
%    1/(mx_max)^2, ...
%    1/(my_max)^2, ...
%   1/(mz_max)^2]);

% Compute the LQR gain
%K_lqr = lqr(A_attitude, B_attitude, Q, R);
%eig(A_attitude - B_attitude * K_lqr);  % Eigenvalue analysis for closed-loop system

%% Longitudinal and Lateral Controller Design (for 100s)
% Longitudinal dynamics
%A_long = A_attitude([2, 4], [2, 4]);
%B_long = B_attitude([2, 4], 2);

% Lateral dynamics
%A_lat = A_attitude([1, 3, 5, 6], [1, 3, 5, 6]);
%B_lat = B_attitude([1, 3, 5, 6], [1, 3]);

% Longitudinal controller
%qmax = deg2rad(1.5);
%alphamax = deg2rad(1);
%my_max = 1600;

%Q_long = diag([1/qmax^2, 1/alphamax^2]);
%R_long = 1/my_max^2;

%K_longitudinal = lqr(A_long, B_long, Q_long, R_long);

% Lateral controller
%pmax = deg2rad(1.5);
%rmax = deg2rad(1.5);
%betamax = deg2rad(1);
%sigmamax = deg2rad(0.5);
%mx_max = 1600;
%mz_max = 1600;

%Q_lat = diag([1/pmax^2, 1/rmax^2, 1/betamax^2, 1/sigmamax^2]);
%R_lat = diag([1/mx_max^2, 1/mz_max^2]);

%K_lateral = lqr(A_lat, B_lat, Q_lat, R_lat);

% References for Lateral and Longitudinal Controllers
%longitudinal_ref = [0; 0];
%lateral_ref = [0; 0; 0; deg2rad(25)];

% === Eigenvalue Analysis ===
%[V, D] = eig(A_attitude);
%eig_vals = diag(D);
%state_names = {'p', 'q', 'r', 'alpha', 'beta', 'sigma'};

% Eigenvalue reporting (optional)
% for i = 1:length(eig_vals)
%     fprintf('\nEigenvalue %d: %.6f %+.6fi\n', i, real(eig_vals(i)), imag(eig_vals(i)));
%     fprintf('State contributions (|eigenvector|):\n');
%     for j = 1:length(state_names)
%         fprintf('  %-6s : %.4f\n', state_names{j}, abs(V(j,i)));
%     end
% end

%% Gain Scheduling for Points in the Trajectory
times = [100, 200, 300, 400, 500];  % [s]
[Mach_vector, K_long_table, K_lat_table] = computeGainSchedule( ...
    times, out, capsule, rho_0, H_s, Mars_radius, stepTime);

% === Save Gain Schedule Data ===
save(fullfile('Controller/data_sets_control', 'GainSchedule.mat'), 'Mach_vector', 'K_long_table', 'K_lat_table');

% === Longitudinal LUT ===
Mach_Klong = Mach_vector;
K_long_data = K_long_table;
save(fullfile('Controller/data_sets_control', 'KLongitudinalLUT.mat'), 'Mach_Klong', 'K_long_data');

% === Lateral LUT ===
K_lat_data = K_lat_table;

%% Signal Editor: Reference Signals for Control Inputs
time = [0, 150, 300, 450];  % Time array for sigma reference
sigma_deg = [100, 70, 120, 120];  % Corresponding sigma values [degrees]

ts_sigma = timeseries(sigma_deg, time);
ts_sigma.Name = 'sigma_ref';

% Create a Simulink dataset for the signal
signalDataset = Simulink.SimulationData.Dataset;
signalDataset = signalDataset.addElement(ts_sigma);

% Save the dataset
save(fullfile('Controller/data_sets_control', 'sigma_ref_dataset.mat'), 'signalDataset');

%% Intermediate Gain Data for Sigma LUT
Mach_Ksigma = [34.0033; 25; 20; 5; 0];
K_sigmaInt_data = linspace(0, 500000, length(Mach_Ksigma))';

save(fullfile('Controller/data_sets_control', 'KsigmaIntLUT.mat'), 'Mach_Ksigma', 'K_sigmaInt_data');

% Optional: Save Beta LUT (commented out)
% Mach_Kbeta = [34.0033; 25; 20; 5; 0];
% K_betaInt_data = linspace(0, 1000, length(Mach_Ksigma))';
% save(fullfile('data_sets_control', 'KbetaIntLUT.mat'), 'Mach_Kbeta', 'K_betaInt_data');

