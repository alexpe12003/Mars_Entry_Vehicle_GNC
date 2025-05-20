clc;

state = out.state.signals.values;
D_L_S = out.D_L_S.signals.values;
F_g = out.F_g.signals.values;
mach = out.mach.signals.values;
p_dyna = out.q.signals.values;

time_step = 200; % [s]
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

A_attitude = A(4:9,4:9)
B_attitude = B(4:9,:)

% Eigenvalues of A
disp('--- Eigenvalues of A ---');
eig_A = eig(A_attitude)

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

K_lqr=lqr(A_attitude,B_attitude,Q,R)
eig(A_attitude- B_attitude*K_lqr)

