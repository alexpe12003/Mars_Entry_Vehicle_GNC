function gain_table = compute_F1_F2_F3_from_costates( ...
    v_star, gamma_star, h_star, hdot_star, D_star, Mars_radius, m, ...
    Cd, Cl, rho, gs, Sref, sigma, L_star)

% === Constants ===
h_s = 7050;  % Mars scale height [m]

% === Time base ===
N = length(v_star);
t = linspace(0, 1, N);
t_rev = fliplr(t);

% === Flip all vectors ===
v_rev     = flipud(v_star);
gamma_rev = flipud(gamma_star);
h_rev     = flipud(h_star);
D_rev     = flipud(D_star);
rho_rev   = flipud(rho);
gs_rev    = flipud(gs);
sigma_rev = flipud(sigma);
L_rev     = flipud(L_star);
Cd_rev = flipud(Cd);
Cl_rev = flipud(Cl);

% === Interpolants ===
f_v     = @(tq) interp1(t_rev, v_rev, tq, 'linear', 'extrap');
f_h     = @(tq) interp1(t_rev, h_rev, tq, 'linear', 'extrap');
f_gam   = @(tq) interp1(t_rev, gamma_rev, tq, 'linear', 'extrap');
f_D     = @(tq) interp1(t_rev, D_rev, tq, 'linear', 'extrap');
f_rho   = @(tq) interp1(t_rev, rho_rev, tq, 'linear', 'extrap');
f_gs    = @(tq) interp1(t_rev, gs_rev, tq, 'linear', 'extrap');
f_sigma = @(tq) interp1(t_rev, sigma_rev, tq, 'linear', 'extrap');
f_L     = @(tq) interp1(t_rev, L_rev, tq, 'linear', 'extrap');
f_Cd     = @(tq) interp1(t_rev, Cd_rev, tq, 'linear', 'extrap');
f_Cl     = @(tq) interp1(t_rev, Cl_rev, tq, 'linear', 'extrap');

% === Adjoint system ===
function dL = adjoint_odes(t, lambda)
    lambda_s  = lambda(1);
    lambda_v  = lambda(2);
    lambda_g  = lambda(3);
    lambda_h  = lambda(4);
    lambda_u  = lambda(5);

    v     = f_v(t);
    h     = f_h(t);
    gamma = f_gam(t);
    D     = f_D(t);
    rho   = f_rho(t);
    g_s   = f_gs(t);
    sigma = f_sigma(t);
    L     = f_L(t);
    cl = f_Cl(t);
    cd = f_Cd(t);

    re_h = Mars_radius + h;
    cosg = cos(gamma);
    sing = sin(gamma);
    cossigma = cos(sigma);

    % Equations from paper (corrected)
    dlambda_s = 0;
    
    
    
    lamGAMdot = -g*lamGAM*sin(gam)/v + g*lamV*cos(gam) + lamGAM*v*sin(gam)/r - lamH*v*cos(gam) + lamS*v*sin(gam)
    lamUdot = LD*lamGAM*rho*v*sin(u)/(2*beta)

    dlambda_v = -cosg * lambda_s ...
        + (rho * v * cd * Sref / m) * lambda_v ...
        + ((rho * cl * cossigma * Sref / (2 * m)) + (cosg / re_h) + (g_s * cosg / v^2)) * lambda_g ...
        - sing * lambda_h;
    
    lamVdot = D_m*LD*lamGAM*cos(u)/v**2 - LD*lamGAM*rho*cos(u)/beta - g*lamGAM*cos(gam)/v**2 - lamGAM*cos(gam)/r - lamH*sin(gam) - lamS*cos(gam) + lamV*rho*v/beta

    dlambda_g = v * sing * lambda_s ...
        + g_s * cosg * lambda_g ...
        + ((v / re_h - g_s / v) * sing) * lambda_g ...
        - v * cosg * lambda_h;

    dlambda_h = (-D / (m * h_s)) * lambda_v ...
        + (L * cossigma / (m * h_s * v) + v * cosg / re_h^2) * lambda_g;
    
    lamHdot = D_m*LD*lamGAM*cos(u)/(H*v) - D_m*lamV/H + lamGAM*v*cos(gam)/r**2
    D_m = rho * V2 / (2 * beta)

    dlambda_u = (-D / (m * v)) * lambda_v;

    dL = [dlambda_s; dlambda_v; dlambda_g; dlambda_h; dlambda_u];
end

% === Boundary conditions ===
gamma_f = gamma_rev(1);
lambda_final = [1; 0; 0; -cot(gamma_f); 0];

% === Integrate backward in time ===
[~, Lambda_rev] = ode45(@adjoint_odes, t_rev, lambda_final);
Lambda = flipud(Lambda_rev);

% === Extract costates ===
lambda_s  = Lambda(:,1);
lambda_v  = Lambda(:,2);
lambda_g  = Lambda(:,3);
lambda_h  = Lambda(:,4);
lambda_u  = Lambda(:,5);

% === Gains ===
F1 = -m .* h_s ./ D_star .* lambda_h;
F2 = lambda_g ./ (v_star .* cos(gamma_star));
F3 = lambda_u;

% === Output ===
gain_table.v = v_star;
gain_table.F1 = F1;
gain_table.F2 = F2;
gain_table.F3 = F3;

save('apollo_gains_final_lift_input.mat', 'gain_table');
end


