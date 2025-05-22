function [Mach_vector, K_long_table, K_lat_table] = computeGainSchedule(times, out, capsule, rho_0, H_s, Mars_radius,stepTime)

% Preallocate outputs
K_long_table = zeros(length(times), 2);       % [N x 2] for q and alpha → My
K_lat_table  = zeros(length(times), 2, 4);    % [N x 2 x 4] for p, r, beta, sigma → Mx, Mz
Mach_vector  = zeros(length(times), 1);       % [N x 1] for interpolation

% Loop through each time to compute LQR gains
for i = 1:length(times)
    t = times(i);
    idx = round(t / stepTime);

    % Extract data
    state     = out.state.signals.values;
    D_L_S     = out.D_L_S.signals.values;
    F_g       = out.F_g.signals.values;
    mach      = out.mach.signals.values;
    p_dyna    = out.q.signals.values;

    % Store Mach number at this time
    Mach_vector(i) = mach(idx);

    % Linearize system
    [A, B, ~, ~] = capsule_linearize(capsule, idx, state, D_L_S, F_g, mach, p_dyna, rho_0, H_s, Mars_radius);
    A_att = A(4:9,4:9);     % attitude states: p, q, r, alpha, beta, sigma
    B_att = B(4:9,:);       % Mx, My, Mz

    %% Longitudinal subsystem (q and alpha → My)
    A_long = A_att([2, 4], [2, 4]);
    B_long = B_att([2, 4], 2);   % My only

    % Tuning parameters
    qmax     = deg2rad(1.5);    % rad/s
    alphamax = deg2rad(1);      % rad
    my_max   = 1600;            % Nm

    Q_long = diag([1/qmax^2, 1/alphamax^2]);
    R_long = 1 / my_max^2;

    % Compute LQR gain
    K_long = lqr(A_long, B_long, Q_long, R_long);
    K_long_table(i, :) = K_long;

    %% Lateral subsystem (p, r, beta, sigma → Mx, Mz)
    A_lat = A_att([1, 3, 5, 6], [1, 3, 5, 6]);
    B_lat = B_att([1, 3, 5, 6], [1, 3]);   % Mx and Mz only

    % Tuning parameters
    pmax      = deg2rad(1.5);    % rad/s
    rmax      = deg2rad(1.5);    % rad/s
    betamax   = deg2rad(1);      % rad
    sigmamax  = deg2rad(0.5);    % rad
    mx_max    = 1600;            % Nm
    mz_max    = 1600;            % Nm

    Q_lat = diag([1/pmax^2, 1/rmax^2, 1/betamax^2, 1/sigmamax^2]);
    R_lat = diag([1/mx_max^2, 1/mz_max^2]);

    % Compute LQR gain
    K_lat = lqr(A_lat, B_lat, Q_lat, R_lat);
    K_lat_table(i,:,:) = K_lat;
    
end
save('GainSchedule.mat', 'Mach_vector', 'K_long_table', 'K_lat_table');
end

