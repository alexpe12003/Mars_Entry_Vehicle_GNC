clc;

load("Trajectory_openLoop_cleaned.mat")

state = out.state.signals.values;
D_L_S = out.D_L_S.signals.values;
F_g = out.F_g.signals.values;
mach = out.mach.signals.values;
p_dyna = out.q.signals.values;

time_step = 100; % [s]
idx = time_step/stepTime;

capsule.m0 = m_0;
capsule.Ixx = I_xx;
capsule.Iyy = I_yy;
capsule.Izz = I_zz;
capsule.Ixy = I_xy;
capsule.Iyz = I_yz;
capsule.Ixz = I_xz;
capsule.S_ref = S_ref;
capsule.d_ref = d_ref;

[A, B, ~, ~] = capsule_linearize(capsule, idx, state, D_L_S, F_g, mach, p_dyna, rho_0, H_s, Mars_radius);

A_attitude = A(4:9,4:9);
B_attitude = B(4:9,:);

% Eigenvalues of A
disp('--- Eigenvalues of A ---');
eig_A = eig(A_attitude);

% Controllability matrix and rank
disp('--- Controllability Analysis ---');
Co = ctrb(A_attitude, B_attitude);
rank_Co = rank(Co);
disp(['Rank of controllability matrix: ', num2str(rank_Co)]);

% State constraints
pmax = deg2rad(3);    
qmax = deg2rad(3);    
rmax = deg2rad(2);    
alfamax = deg2rad(3); 
betamax = deg2rad(3); 
sigmamax = deg2rad(3);

% Actuator limits (moments)
mx_max = 1000; 
my_max = 1000; 
mz_max = 1000; 

% Build Q and R matrices
Q = diag([ ...
    1/(pmax)^2, ...
    1/(qmax)^2, ...
    1/(rmax)^2, ...
    1/(alfamax)^2, ...
    1/(betamax)^2, ...
    1/(sigmamax)^2]);

R = diag([ ...
    1/(mx_max)^2, ...
    1/(my_max)^2, ...
    1/(mz_max)^2]);

K_lqr = lqr(A_attitude,B_attitude,Q,R);
eig(A_attitude - B_attitude*K_lqr);

%% Lateral and Longitudinal Controller

A_long = A_attitude([2, 4], [2, 4]);
B_long = B_attitude([2, 4], 2);

A_lat = A_attitude([1, 3, 5, 6], [1, 3, 5, 6]);
B_lat = B_attitude([1, 3, 5, 6], [1, 3]);

qmax = deg2rad(1.5);
alphamax = deg2rad(1);
my_max = 1600;

Q_long = diag([1/qmax^2, 1/alphamax^2]);
R_long = 1/my_max^2;

K_longitudinal = lqr(A_long, B_long, Q_long, R_long);

pmax = deg2rad(1.5);
rmax = deg2rad(1.5);
betamax = deg2rad(1);
sigmamax = deg2rad(0.5);
mx_max = 1600;
mz_max = 1600;

Q_lat = diag([1/pmax^2, 1/rmax^2, 1/betamax^2, 1/sigmamax^2]);
R_lat = diag([1/mx_max^2, 1/mz_max^2]);

K_lateral = lqr(A_lat, B_lat, Q_lat, R_lat);

longitudinal_ref = [0; 0];
lateral_ref = [0; 0; 0; deg2rad(25)];

% Eigenvalue analysis
[V, D] = eig(A_attitude);
eig_vals = diag(D);
state_names = {'p', 'q', 'r', 'alpha', 'beta', 'sigma'};

fprintf('\n======= Eigenvalue Analysis =======\n');
for i = 1:length(eig_vals)
    fprintf('\nEigenvalue %d: %.6f %+.6fi\n', i, real(eig_vals(i)), imag(eig_vals(i)));
    fprintf('State contributions (|eigenvector|):\n');
    for j = 1:length(state_names)
        fprintf('  %-6s : %.4f\n', state_names{j}, abs(V(j,i)));
    end
end

%% Gain Schedule
times = [100 200 300 400 500];  % [s]
[Mach_vector, K_long_table, K_lat_table] = computeGainSchedule( ...
    times, out, capsule, rho_0, H_s, Mars_radius, stepTime);

% Ensure data_sets directory exists
if ~exist('data_sets', 'dir')
    mkdir('data_sets');
end

% Save gain schedule data
save(fullfile('Controller/data_sets', 'GainSchedule.mat'), 'Mach_vector', 'K_long_table', 'K_lat_table');

Mach_Klong = Mach_vector;
K_long_data = K_long_table;
save(fullfile('Controller/data_sets', 'KLongitudinalLUT.mat'), 'Mach_Klong', 'K_long_data');

K_lat_data = K_lat_table;

%% Signal Editor
time = [0, 150, 300, 450];
sigma_deg = [100, 70, 120, 120];

ts_sigma = timeseries(sigma_deg, time);
ts_sigma.Name = 'sigma_ref';

signalDataset = Simulink.SimulationData.Dataset;
signalDataset = signalDataset.addElement(ts_sigma);

save(fullfile('Controller/data_sets', 'sigma_ref_dataset.mat'), 'signalDataset');

Mach_Ksigma = [34.0033; 25; 20; 5; 0];
K_sigmaInt_data = linspace(0, 500000, length(Mach_Ksigma))';
Mach_Kbeta = [34.0033; 25; 20; 5; 0];
K_betaInt_data = linspace(0, 1000, length(Mach_Ksigma))';

save(fullfile('Controller/data_sets', 'KsigmaIntLUT.mat'), 'Mach_Ksigma', 'K_sigmaInt_data');
save(fullfile('Controller/data_sets', 'KbetaIntLUT.mat'), 'Mach_Kbeta', 'K_betaInt_data');
