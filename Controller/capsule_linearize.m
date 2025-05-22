function [A, B, state_0, state_dot_0] = capsule_linearize(capsule, idx, state, D_L_S, F_g, mach, p_dyna, rho_0, H_s, Mars_radius)

% Vehicle properties
    m = capsule.m0; % [kg]
    I_xx = capsule.Ixx; % [kg/m^2]
    I_yy = capsule.Iyy; % [kg/m^2]
    I_zz = capsule.Izz; % [kg/m^2]
    I_xy = capsule.Ixy; % [kg/m^2]
    I_yz = capsule.Iyz; % [kg/m^2]
    I_xz = capsule.Ixz; % [kg/m^2]
    Sref = capsule.S_ref;
    dref = capsule.d_ref;

    % Denominator term
    I_star = I_xx * I_zz - I_xz^2;
    % Coefficients used in linearized dynamics
    Iq1 = I_xz / I_yy;
    Iq2 = (I_zz - I_xx) / I_yy;
    Ir1 = ((I_xx - I_yy) * I_xx - I_xz^2) / I_star;
    Ir2 = ((-I_xx + I_yy - I_zz) * I_xz) / I_star;
    Ip1 = ((I_xx - I_yy + I_zz) * I_xz) / I_star;
    Ip2 = ((I_yy - I_zz) * I_zz - I_xz^2) / I_star;
    gamma_dot_0 = 1;

% Assign variables to state variables
    V0 = state(idx, 7);
    gamma0 = state(idx, 8);
    R0 = state(idx, 1);
    alpha0 = state(idx, 4);
    alpha0 = deg2rad(-23.82);
    sigma0 = state(idx, 6);

    D0 = D_L_S(idx, 1);
    L0 = D_L_S(idx, 2);

    p_dyna0 = p_dyna(idx);

    M0 = mach(idx);
    g0 = F_g(idx)/m;

% Compute equilibrium angular rates
    c1 = g0/V0 * cos(gamma0) * sin(sigma0);
    c2 = L0/m/V0 * tan(gamma0) * sin(sigma0);
    p0 = c1*sin(alpha0) + c2*cos(alpha0);
    q0 = L0/m/V0 - g0/V0 * cos(gamma0) * cos(sigma0);
    r0 = -c1 * cos(alpha0) + c2 * sin(alpha0);

% Compute local aerodynamic derivatives
    [alphaT_deg, CD0, CL0, CS0, Cm0, Cn0, Cmq0, dCm_dalpha, dCL_dalpha, dCD_dalpha, dCL_dbeta, dCS_dbeta, dCn_dbeta, dCD_dM, dCl_dbeta, dCL_dM, dCm_dM] = computeAerodynamicDerivatives(M0, alpha0, 0, rho_0, V0);

    drho_dh = -rho_0/H_s * exp(-(R0-Mars_radius)/H_s);
    da_dh = 0;
%% A matrix
% Row 1: ΔV̇
    a_VV      = -(1 / (m * V0)) * (M0 * dCD_dM * p_dyna0 * Sref + 2 * D0);
    a_Vgamma  = -g0 * cos(gamma0);
    a_VR      = (2 * g0 / R0) * sin(gamma0) ...
              - (1 / (m * V0)) * (dCD_dM * M0^2 * da_dh * p_dyna0 * Sref) ...
              + (0.5 / m) * CD0 * Sref * V0^2 * drho_dh;
    a_Valpha  = -(1 / m) * dCD_dalpha * p_dyna0 * Sref;
    
% Row 2: Δγ̇
    a_gammaV      = (1 / V0) * (-gamma_dot_0 + (2 * V0 / R0) * cos(gamma0)) ...
                  + (cos(sigma0) / (m * V0^2)) * (M0 * dCL_dM * p_dyna0 * Sref + 2 * L0);
    a_gammagamma  = - (V0 / R0 - g0 / V0) * sin(gamma0);
    a_gammaR      = ((2 * g0 / V0) - (V0 / R0)) * cos(gamma0) / R0 ...
                  + (cos(sigma0) / (m * V0)) * ...
                    (0.5 * CL0 * Sref * (V0^2) * drho_dh ...
                    - dCL_dM * (M0^2/V0) * da_dh * p_dyna0 * Sref);
    a_gammaalpha  = (cos(sigma0) / (m * V0)) * dCL_dalpha * p_dyna0 * Sref;
    a_gammabeta   = -(sin(sigma0) / (m * V0)) * dCS_dbeta * p_dyna0 * Sref;
    a_gammasigma  = - (L0 / (m * V0)) * sin(sigma0);
    
% Row 3: ΔṘ
    a_RV     = sin(gamma0);
    a_Rgamma = V0 * cos(gamma0);
    
% Row 4: Δṗ
    a_pp     = Ip1 * q0;
    a_pq     = Ip1 * p0 + Ip2 * r0;
    a_pr     = Ip2 * q0;
    a_pbeta  = (I_zz / I_star) * dCl_dbeta * p_dyna0 * Sref * dref ...
             + (I_xz / I_star) * dCn_dbeta * p_dyna0 * Sref * dref;
    
% Row 5: Δq̇
    a_qV      = (1 / I_yy) * M0 * V0 * dCm_dM * p_dyna0 * Sref * dref;
    a_qR      = -(1 / I_yy) * M0^2 * V0 * dCm_dM * p_dyna0 * Sref * dref * da_dh;
    a_qp      = -2 * Iq1 * p0 + Iq2 * r0;
    a_qq      = 0;  % If no pitch damping coefficient available
    a_qr      = Iq2 * p0 + 2 * Iq1 * r0;
    a_qalpha  = (1 / I_yy) * dCm_dalpha * p_dyna0 * Sref * dref;
    a_qbeta   = 0;  % Set if no beta moment derivative known
    
% Row 6: Δṙ
    a_rp     = Ir1 * p0;
    a_rq     = Ir1 * p0 + Ir2 * r0;
    a_rr     = Ir2 * q0;
    a_rbeta  = (I_xx / I_star) * dCn_dbeta * p_dyna0 * Sref * dref ...
             + (I_xz / I_star) * dCl_dbeta * p_dyna0 * Sref * dref;
    
%Row 7: Δalfa
    a_alphaV = - (g0 / V0^2) * cos(gamma0) * cos(sigma0) ...
             - (1 / (m * V0^2)) * (M0 * dCL_dM + CL0) * p_dyna0 * Sref;
    a_alphagamma = - (g0 / V0) * sin(gamma0) * cos(sigma0);
    a_alphaR = - (2 * g0 / (R0 * V0)) * cos(gamma0) * cos(sigma0) ...
             + (1 / (m * V0)) * ...
               (0.5 * CL0 * Sref * (V0^2) * drho_dh ...
              - dCL_dM * (M0^2/V0) * da_dh * p_dyna0 * Sref);
    a_alphaq = 1;
    a_alphaalpha = - (1 / (m * V0)) * dCL_dalpha * p_dyna0 * Sref;
    a_alphasigma = - (g0 / V0) * cos(gamma0) * sin(sigma0);
    
%Row 8: delta beta
    a_betaV = (g0 / V0^2) * cos(gamma0) * sin(sigma0);
    a_betagamma = (g0 / V0) * sin(gamma0) * sin(sigma0);
    a_betaR = (2 * g0 / (R0 * V0)) * cos(gamma0) * sin(sigma0);
    a_betap = sin(alpha0);
    a_betar = -cos(alpha0);
    a_betabeta = - (1 / (m * V0)) * dCS_dbeta * p_dyna0 * Sref;
    a_betasigma = - (g0 / V0) * cos(gamma0) * cos(sigma0);
    
%Row 9: delta sigma
    a_sigmaV = (tan(gamma0) * sin(sigma0) / (m * V0^2)) ...
             * (M0 * dCL_dM + CL0) * p_dyna0 * Sref;
    a_sigmagamma = (L0 / (m * V0)) * sin(sigma0);
    a_sigmaR = (tan(gamma0) * sin(sigma0) / (m * V0)) * ...
              (0.5 * CL0 * Sref * (V0^2) * drho_dh ...
             - dCL_dM * (M0^2/V0) * da_dh * p_dyna0 * Sref);
    a_sigmap = -cos(alpha0);
    a_sigmar = -sin(alpha0);
    a_sigmaalpha = (tan(gamma0) * sin(sigma0) / (m * V0)) * dCL_dalpha * p_dyna0 * Sref;
    a_sigmabeta = (tan(gamma0) * cos(sigma0) / (m * V0)) * dCS_dbeta * p_dyna0 * Sref ...
                - (L0 / (m * V0)) ...
                + (g0 / V0) * cos(gamma0) * cos(sigma0);
    a_sigmasigma = (tan(gamma0) * cos(sigma0) / (m * V0)) * L0;

%% B matrix
    b_px = I_zz / I_star;
    b_pz = I_xz / I_star;
    b_qy = 1 / I_yy;
    b_rx = I_xz / I_star;
    b_rz = I_xx / I_star;

%% === SYSTEM MATRICES ===

% Initialize A Matrix (9x9)
A = zeros(9,9);

% Row 1: ΔV̇
A(1,1) = a_VV;       A(1,2) = a_Vgamma;     A(1,3) = a_VR;    A(1,7) = a_Valpha;

% Row 2: Δγ̇
A(2,1) = a_gammaV;   A(2,2) = a_gammagamma; A(2,3) = a_gammaR;
A(2,7) = a_gammaalpha; A(2,8) = a_gammabeta; A(2,9) = a_gammasigma;

% Row 3: ΔṘ
A(3,1) = a_RV;       A(3,2) = a_Rgamma;

% Row 4: Δṗ
A(4,4) = a_pp;       A(4,5) = a_pq;         A(4,6) = a_pr;        A(4,8) = a_pbeta;

% Row 5: Δq̇
A(5,1) = a_qV;       A(5,3) = a_qR;         A(5,4) = a_qp;
A(5,5) = a_qq;       A(5,6) = a_qr;         A(5,7) = a_qalpha;    A(5,8) = a_qbeta;

% Row 6: Δṙ
A(6,4) = a_rp;       A(6,5) = a_rq;         A(6,6) = a_rr;        A(6,8) = a_rbeta;

% Remaining rows (alpha, beta, sigma) would follow here...
% e.g., A(7,:), A(8,:), A(9,:) with their respective derivatives

%Row 7: Δalfa
A(7,1) = a_alphaV; A(7,2) = a_alphagamma; A(7,3) = a_alphaR;
A(7,5) = a_alphaq; A(7,7) = a_alphaalpha; A(7,9) = a_alphasigma;

%Row 8: delta beta
A(8,1) = a_betaV; A(8,2) = a_betagamma; A(8,3) = a_betaR;
A(8,4) = a_betap;
A(8,6) = a_betar;
A(8,8) = a_betabeta;
A(8,9) = a_betasigma;

%Row 9: delta sigma
A(9,1) = a_sigmaV;
A(9,2) = a_sigmagamma;
A(9,3) = a_sigmaR;
A(9,4) = a_sigmap;
A(9,6) = a_sigmar;
A(9,7) = a_sigmaalpha;
A(9,8) = a_sigmabeta;
A(9,9) = a_sigmasigma;

% Initialize B Matrix (9x3)
B = zeros(9,3);
B(4,1) = b_px;
B(4,3) = b_pz;
B(5,2) = b_qy;
B(6,1) = b_rx;
B(6,3) = b_rz;

state_0 = 69;
state_dot_0 = 69;

end