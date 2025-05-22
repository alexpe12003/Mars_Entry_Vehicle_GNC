clc;

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

% State constraints (define in rad or proper units)
%vmax = 100;             % [m/s] (assuming V)
%gammamax = deg2rad(3);% [rad]
%Rmax = deg2rad(4);    % [rad]
pmax = deg2rad(3);    % [rad/s]
qmax = deg2rad(3);    % [rad/s]
rmax = deg2rad(2);    % [rad/s]
alfamax = deg2rad(3); % [rad]
betamax = deg2rad(3); % [rad]
sigmamax = deg2rad(3);% [rad]

% Actuator limits (moments)
mx_max = 1000; % [Nm]
my_max = 1000; % [Nm]
mz_max = 1000; % [Nm]


% Build Q (penalize deviation relative to max allowed)
Q = diag([ ...
    1/(pmax)^2, ...
    1/(qmax)^2, ...
    1/(rmax)^2, ...
    1/(alfamax)^2, ...
    1/(betamax)^2, ...
    1/(sigmamax)^2]);

% Build R (penalize control effort relative to actuator max)
R = diag([ ...
    1/(mx_max)^2, ...
    1/(my_max)^2, ...
    1/(mz_max)^2]);

K_lqr=lqr(A_attitude,B_attitude,Q,R);
eig(A_attitude- B_attitude*K_lqr);

%% Lateral and Longitudinal Controller

% Indices: [p=1, q=2, r=3, alpha=4, beta=5, sigma=6]
% Longitudinal subsystem (q and alpha)
A_long = A_attitude([2, 4], [2, 4]);
B_long = B_attitude([2, 4], 2);   % My only

% Lateral subsystem (p, r, beta, sigma)
A_lat = A_attitude([1, 3, 5, 6], [1, 3, 5, 6]);
B_lat = B_attitude([1, 3, 5, 6], [1, 3]);   % Mx and Mz only

qmax = deg2rad(1.5);       % rad/s
alphamax = deg2rad(1);   % rad
my_max = 1600;           % Nm

Q_long = diag([1/qmax^2, 1/alphamax^2]);
R_long = 1/my_max^2;

K_longitudinal = lqr(A_long, B_long, Q_long, R_long);

pmax = deg2rad(1.5);     % rad/s
rmax = deg2rad(1.5);     % rad/s
betamax = deg2rad(1);  % rad
sigmamax = deg2rad(1); % rad
mx_max = 1600;
mz_max = 1600;

Q_lat = diag([1/pmax^2, 1/rmax^2, 1/betamax^2, 1/sigmamax^2]);
R_lat = diag([1/mx_max^2, 1/mz_max^2]);

K_lateral = lqr(A_lat, B_lat, Q_lat, R_lat);

% === 1. Yaw Controller: [r, beta] → ΔMz ===
A_rbeta = A_attitude([3, 5], [3, 5]);
B_rbeta = B_attitude([3, 5], 3);  % Mz

rmax     = deg2rad(1.5);   % rad/s
betamax  = deg2rad(1);     % rad
mz_max   = 1600;           % Nm

Q_rbeta = diag([1/rmax^2, 1/betamax^2]);
R_rbeta = 1/mz_max^2;

K_rbeta = lqr(A_rbeta, B_rbeta, Q_rbeta, R_rbeta);

% === 2. Roll-Bank Controller: [p, sigma] → ΔMx ===
A_psigma = A_attitude([1, 6], [1, 6]);
B_psigma = B_attitude([1, 6], 1);  % Mx

pmax      = deg2rad(1.5);   % rad/s
sigmamax  = deg2rad(4);     % rad
mx_max    = 1600;           % Nm

Q_psigma = diag([1/pmax^2, 1/sigmamax^2]);
R_psigma = 1/mx_max^2;

K_psigma = lqr(A_psigma, B_psigma, Q_psigma, R_psigma);

longitudinal_ref = [0;0];
lateral_ref = [0;0;0;deg2rad(25)];

% Compute eigenvalues and eigenvectors of A_attitude
[V, D] = eig(A_attitude);

% Extract eigenvalues (on diagonal of D)
eig_vals = diag(D);

% Define state names in attitude block
state_names = {'p', 'q', 'r', 'alpha', 'beta', 'sigma'};

fprintf('\n======= Eigenvalue Analysis =======\n');
for i = 1:length(eig_vals)
    fprintf('\nEigenvalue %d: %.6f %+.6fi\n', i, real(eig_vals(i)), imag(eig_vals(i)));
    fprintf('State contributions (|eigenvector|):\n');
    
    for j = 1:length(state_names)
        fprintf('  %-6s : %.4f\n', state_names{j}, abs(V(j,i)));
    end
end
