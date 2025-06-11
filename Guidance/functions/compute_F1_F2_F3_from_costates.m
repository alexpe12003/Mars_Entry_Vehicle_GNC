function gain_table = compute_F1_F2_F3_from_costates( ...
    v_star, gamma_star, h_star, hdot_star, D_star, ...
    Mars_radius, m, rho, gs, Sref, sigma, ballistic_beta, L_over_D,t)

% Inputs:
% v_star, gamma_star, h_star, hdot_star, D_star - State vectors for velocity, flight path angle, altitude, altitude rate, and drag
% Mars_radius - Radius of Mars [m]
% m - Mass of the vehicle [kg]
% rho - Atmospheric density [kg/m^3]
% gs - Gravitational acceleration per unit mass [m/s^2]
% Sref - Reference area [m^2]
% sigma - Bank angle [rad]
% ballistic_beta - Ballistic coefficient [kg/m^2]
% L_over_D - Lift-to-drag ratio [dimensionless]

% === Constants ===
h_s = 7050; % A scaling factor for the altitude
N = length(v_star); % Number of data points

% === Time vector and reverse time for integration ===
t = linspace(0, 1, N);
t_rev = fliplr(t);  % Reversed time vector

% === Reverse the state vectors ===
v_rev     = flipud(v_star);
gamma_rev = flipud(gamma_star);
h_rev     = flipud(h_star);
D_rev     = flipud(D_star);
rho_rev   = flipud(rho);
gs_rev    = flipud(gs);
sigma_rev = flipud(sigma);

% === Create interpolant functions for state variables ===
% These functions will be used to interpolate the state vectors at any time point
f_v     = @(tq) interp1(t_rev, v_rev, tq, 'linear', 'extrap');
f_h     = @(tq) interp1(t_rev, h_rev, tq, 'linear', 'extrap');
f_gam   = @(tq) interp1(t_rev, gamma_rev, tq, 'linear', 'extrap');
f_D     = @(tq) interp1(t_rev, D_rev, tq, 'linear', 'extrap');
f_rho   = @(tq) interp1(t_rev, rho_rev, tq, 'linear', 'extrap');
f_gs    = @(tq) interp1(t_rev, gs_rev, tq, 'linear', 'extrap');
f_sigma = @(tq) interp1(t_rev, sigma_rev, tq, 'linear', 'extrap');

% === Adjoint ODEs (Costates) ===
% Define the system of adjoint equations
function dL = adjoint_odes(t, lambda)
    % Unpack lambda (costates)
    lambda_s  = lambda(1);
    lambda_v  = lambda(2);
    lambda_g  = lambda(3);
    lambda_h  = lambda(4);
    lambda_u  = lambda(5);

    % Interpolate the state variables at the current time
    v     = f_v(t);
    h     = f_h(t);
    gamma = f_gam(t);
    D     = f_D(t);
    rho   = f_rho(t);
    g_s   = f_gs(t);
    sigma = f_sigma(t);

    % Compute intermediate variables for the adjoint system
    re_h = Mars_radius + h;
    cosg = cos(gamma);
    sing = sin(gamma);
    cossigma = cos(sigma);

    % === Compute Cd, Cl, and L from D, rho, beta, and L/D ===
    Cd = m / (ballistic_beta * Sref);  % Drag coefficient
    L  = L_over_D * D;                 % Lift
    Cl = 2 * L / (rho * v^2 * Sref);   % Lift coefficient

    % Adjoint system equations (costate equations)
    dlambda_s = 0;

    dlambda_v = -cosg * lambda_s ...
        + (rho * v * Cd * Sref / m) * lambda_v ...
        + ((rho * Cl * cossigma * Sref / (2 * m)) + (cosg / re_h) + (g_s * cosg / v^2)) * lambda_g ...
        - sing * lambda_h;

    dlambda_g = v * sing * lambda_s ...
        + g_s * cosg * lambda_v ...
        + ((v / re_h - g_s / v) * sing) * lambda_g ...
        - v * cosg * lambda_h;

    dlambda_h = (-D / (m * h_s)) * lambda_v ...
        + (L * cossigma / (m * h_s * v) + v * cosg / re_h^2) * lambda_g;

    dlambda_u = (-D / (m * v)) * lambda_g;

    % Return the derivatives of the costates
    dL = [dlambda_s; dlambda_v; dlambda_g; dlambda_h; dlambda_u];
end

% === Final costate values and integration ===
gamma_f = gamma_rev(1);  % Final flight path angle
lambda_final = [1; 0; 0; -cot(gamma_f); 0];  % Initial final costates

% Integrate the adjoint ODEs using ode45
[~, Lambda_rev] = ode45(@adjoint_odes, t_rev, lambda_final);

% Reverse the solution to align with the original time vector
Lambda = flipud(Lambda_rev);

% Extract the individual costates
lambda_h = Lambda(:, 4);
lambda_g = Lambda(:, 3);
lambda_u = Lambda(:, 5);

% === Compute F1, F2, F3 ===
F1 = -m .* h_s .* lambda_h ./ D_star;  % F1 - Scaling based on drag force
F2 = lambda_g ./ (v_star .* cos(gamma_star));  % F2 - Based on velocity and flight path angle
F3 = lambda_u;  % F3 - Based on the adjoint variable for control input

% === Create the gain table ===
gain_table.v  = v_star;
gain_table.F1 = F1;
gain_table.F2 = F2;
gain_table.F3 = F3;

end




