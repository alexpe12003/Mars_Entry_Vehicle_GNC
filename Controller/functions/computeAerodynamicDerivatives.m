function [alphaT_deg, C_D, C_L, C_S, C_m, C_n, Cm_q, dCm_dalpha, dCL_dalpha, dCD_dalpha, dCL_dbeta, dCS_dbeta, dCn_dbeta, dCD_dM, dCl_dbeta, dCL_dM, dCm_dM] = computeAerodynamicDerivatives(M0, alpha_rad, beta_rad, rho0, V0)
% Computes all Apollo aerodynamic coefficients and their derivatives
% Inputs in degrees, outputs in coefficient units and per-radian derivatives
Sref = 12; % [m^2]
r_cg = [-0.137; 0; 1.8];  % [m] in body frame
cref = 3.9; % [m^2] - Vehicle diameter
%% === Convert angles to deegres ===
alpha = rad2deg(alpha_rad);
beta  = rad2deg(beta_rad);
cosA = cosd(alpha); sinA = sind(alpha);
cosB = cosd(beta);  sinB = sind(beta);

%% === Compute total angle of attack alpha_T ===
cosT = cosA * cosB;
alphaT_deg = acosd(cosT);
sinT = sind(alphaT_deg);
epsilon = 1e-6;
if sinT < epsilon
    warning("AlphaT near singularity, using epsilon for stability");
    sinT = epsilon;
end

%% === Load Lookup Tables ===
data = load('ApolloDragCoeff.mat');
CD_table = data.CD_T;
Mach_numbers_CD = data.Mach_numbers_CD;
alpha_T_CD = data.alpha_T_CD;

data = load('ApolloLiftCoeff.mat');
CLT_table = data.CL_T;
Mach_numbers_CL = data.Mach_numbers_CL;
alpha_T_CL = data.alpha_T_CL;

data = load('ApolloPitchMomentCoeff.mat');
CmT_table = data.Cm_T;
Mach_numbers_Cm = data.Mach_numbers_Cm;
alpha_T_Cm = data.alpha_T_Cm;

data = load('ApolloPitchDampingCoeff.mat');
Cm_q_table = data.Cm_qT;
Mach_numbers_Cm_q = data.Mach_numbers_Cm_q;
alpha_T_Cm_q = data.alpha_T_Cm_q;

%% === Build Interpolants ===
%% === Create NDGRID Grids ===
[M_CD, A_CD] = ndgrid(Mach_numbers_CD, alpha_T_CD);
[M_CL, A_CL] = ndgrid(Mach_numbers_CL, alpha_T_CL);
[M_Cm, A_Cm] = ndgrid(Mach_numbers_Cm, alpha_T_Cm);
[M_Cmq, A_Cmq] = ndgrid(Mach_numbers_Cm_q, alpha_T_Cm_q);

%% === Build Interpolants with NDGRID ===
CD_interp   = griddedInterpolant(M_CD, A_CD, CD_table, 'linear', 'nearest');
CLT_interp  = griddedInterpolant(M_CL, A_CL, CLT_table, 'linear', 'nearest');
CmT_interp  = griddedInterpolant(M_Cm, A_Cm, CmT_table, 'linear', 'nearest');
Cm_q_interp = griddedInterpolant(M_Cmq, A_Cmq, Cm_q_table, 'linear', 'nearest');


%% === Interpolated Values ===
CD_T = CD_interp(M0, alphaT_deg);
CL_T = CLT_interp(M0, alphaT_deg);
Cm_T = CmT_interp(M0, alphaT_deg);
Cm_q = Cm_q_interp(M0, alphaT_deg);

%% === Coefficient Transformations ===
C_D = CD_T;
C_L = (sinA / sinT) * CL_T;
C_S = (cosA * sinB / sinT) * CL_T;
C_m = (sinA * cosB / sinT) * Cm_T;
C_n = (-sinB / sinT) * Cm_T;
%% === Derivative d(alphaT)/d(alpha), d(alphaT)/d(beta) ===
dAlphaT_dalpha = -(sinA * cosB) / sinT;
dAlphaT_dbeta  = -(cosA * sinB) / sinT;

%% === Finite differences for CLT, CmT, CDT ===
d_alphaT = 1e-6;
Cm_Tp = CmT_interp(M0, alphaT_deg + d_alphaT);
Cm_Tm = CmT_interp(M0, alphaT_deg - d_alphaT);
dCmT_dalphaT = (Cm_Tp - Cm_Tm) / deg2rad(2*d_alphaT);  % per radian 

CL_Tp = CLT_interp(M0, alphaT_deg + d_alphaT);
CL_Tm = CLT_interp(M0, alphaT_deg - d_alphaT);
dCLT_dalphaT = (CL_Tp - CL_Tm) / deg2rad(2*d_alphaT);  % per radian

CD_Tp = CD_interp(M0, alphaT_deg + d_alphaT);
CD_Tm = CD_interp(M0, alphaT_deg - d_alphaT);
dCDT_dalphaT = (CD_Tp - CD_Tm) / deg2rad(2*d_alphaT);  % per radian

%% Aerodynamic Contribution
% Step size
dalpha = 1e-5;

% Forward/backward alpha
alpha_plus = alpha_rad + dalpha;
alpha_minus = alpha_rad - dalpha;

% Compute alpha_T ± ΔαT in degrees
cosT_p = cos(alpha_plus) * cos(beta_rad);
alphaT_plus = acosd(cosT_p);

cosT_m = cos(alpha_minus) * cos(beta_rad);
alphaT_minus = acosd(cosT_m);

% Recompute force coefficients at α ± Δα
CD_p = CD_interp(M0, alphaT_plus);
CL_p = (sin(alpha_plus) / sin(alphaT_plus)) * CLT_interp(M0, alphaT_plus);
CS_p = (cosT_p * sin(alpha_plus) / sin(alphaT_plus)) * CLT_interp(M0, alphaT_plus);  % If C_S is also a function of α, interpolate it here

% Recompute force coefficients at α ± Δα
CD_m = CD_interp(M0, alphaT_minus);
CL_m = (sin(alpha_minus) / sin(alphaT_minus)) * CLT_interp(M0, alphaT_minus);
CS_m = (cosT_m * sin(alpha_minus) / sin(alphaT_minus)) * CLT_interp(M0, alphaT_minus);

% Convert coefficients to forces
q_dyn = 0.5 * rho0 * V0^2;

D_p = CD_p * q_dyn * Sref;
L_p = CL_p * q_dyn * Sref;
S_p = CS_p * q_dyn * Sref;

D_m = CD_m * q_dyn * Sref;
L_m = CL_m * q_dyn * Sref;
S_m = CS_m * q_dyn * Sref;

% Compute transformed forces at α + Δα
[Xp, Yp, Zp] = aeroFrame2apexFrameForces(D_p, L_p, S_p, alpha_plus, beta_rad);
Fp = [Xp; Yp; Zp];
Mp = cross(r_cg, Fp);
Cm_offset_px = Mp(1) / (q_dyn * Sref * cref);
Cm_offset_py = Mp(2) / (q_dyn * Sref * cref);
Cm_offset_pz = Mp(3) / (q_dyn * Sref * cref);

% Compute transformed forces at α - Δα
[Xm, Ym, Zm] = aeroFrame2apexFrameForces(D_m, L_m, S_m, alpha_minus, beta_rad);
Fm = [Xm; Ym; Zm];
Mm = cross(r_cg, Fm);
Cm_offset_mx = Mm(1) / (q_dyn * Sref * cref);
Cm_offset_my = Mm(2) / (q_dyn * Sref * cref);
Cm_offset_mz = Mm(3) / (q_dyn * Sref * cref);

% Central difference
dCm_offset_dalphax = (Cm_offset_px - Cm_offset_mx) / (2 * dalpha);
dCm_offset_dalphay = (Cm_offset_py - Cm_offset_my) / (2 * dalpha);
dCm_offset_dalphaz = (Cm_offset_pz - Cm_offset_mz) / (2 * dalpha);


%% === dCm/dalpha ===
dCm_dalpha = dCmT_dalphaT * (sinA * cosB / sinT) + Cm_T * (cosA * cosB * sinT - sinA * cosB * cosT) / (sinT^2) + dCm_offset_dalphax + dCm_offset_dalphay + dCm_offset_dalphaz;

%% === dCL/dalpha ===
dCL_dalpha = dCLT_dalphaT * (sinA / sinT) + ...
    CL_T * (cosA * sinT - sinA * cosT) / (sinT^2);

%% === dCD/dalpha ===
dCD_dalpha = dCDT_dalphaT * dAlphaT_dalpha;

%% === dCL/dbeta ===
dCL_dbeta = -CL_T * sinA * cosT * dAlphaT_dbeta / (sinT^2);

%% === dCS/dbeta ===
term1 = (cosA * cosB / sinT) * CL_T;
term2 = -(cosA * sinB * cosT * dAlphaT_dbeta / sinT^2) * CL_T;
dCS_dbeta = term1 + term2;

%% === dCn/dbeta ===
dCn_dbeta = - (cosB / sinT - sinB * cosT * dAlphaT_dbeta / sinT^2)*Cm_T;

% Step size for sideslip
dbeta = 1e-5;

% Forward/backward beta
beta_plus = beta_rad + dbeta;
beta_minus = beta_rad - dbeta;

% Compute cosT for each beta
cosT_bp = cos(alpha_rad) * cos(beta_plus);
cosT_bm = cos(alpha_rad) * cos(beta_minus);
alphaT_plus_beta = acosd(cosT_bp);
alphaT_minus_beta = acosd(cosT_bm);

% Interpolate force coefficients at ±dbeta
CD_bp = CD_interp(M0, alphaT_plus_beta);
CLT_bp = CLT_interp(M0, alphaT_plus_beta);
CL_bp = (sin(alpha_rad) / sind(alphaT_plus_beta)) * CLT_bp;
CS_bp = (cos(alpha_rad) * sin(beta_plus) / sind(alphaT_plus_beta)) * CLT_bp;

CD_bm = CD_interp(M0, alphaT_minus_beta);
CLT_bm = CLT_interp(M0, alphaT_minus_beta);
CL_bm = (sin(alpha_rad) / sind(alphaT_minus_beta)) * CLT_bm;
CS_bm = (cos(alpha_rad) * sin(beta_minus) / sind(alphaT_minus_beta)) * CLT_bm;

% Compute body-frame forces
D_bp = CD_bp * q_dyn * Sref;
L_bp = CL_bp * q_dyn * Sref;
S_bp = CS_bp * q_dyn * Sref;

D_bm = CD_bm * q_dyn * Sref;
L_bm = CL_bm * q_dyn * Sref;
S_bm = CS_bm * q_dyn * Sref;

% Transform to body-frame forces
[Fx_bp, Fy_bp, Fz_bp] = aeroFrame2apexFrameForces(D_bp, L_bp, S_bp, alpha_rad, beta_plus);
Fp_beta = [Fx_bp; Fy_bp; Fz_bp];
Mp_beta = cross(r_cg, Fp_beta);
Cn_offset_p = Mp_beta(3) / (q_dyn * Sref * cref);

[Fx_bm, Fy_bm, Fz_bm] = aeroFrame2apexFrameForces(D_bm, L_bm, S_bm, alpha_rad, beta_minus);
Fm_beta = [Fx_bm; Fy_bm; Fz_bm];
Mm_beta = cross(r_cg, Fm_beta);
Cn_offset_m = Mm_beta(3) / (q_dyn * Sref * cref);

% Central difference
dCn_offset_dbeta = (Cn_offset_p - Cn_offset_m) / (2 * dbeta);

% Add to total
dCn_dbeta = dCn_dbeta + dCn_offset_dbeta;

%% === Finite difference for dCD/dM ===
d_M = 1e-6;  % Small Mach step
CD_TpM = CD_interp(M0 + d_M, alphaT_deg);
CD_TmM = CD_interp(M0 - d_M, alphaT_deg);
dCD_dM = (CD_TpM - CD_TmM) / (2 * d_M); % per Mach

%% === Finite difference for dCl/dbeta ===
dCl_dbeta = 0;

%% === Finite difference for dCL/dM ===
CL_TpM = CLT_interp(M0 + d_M, alphaT_deg);
CL_TmM = CLT_interp(M0 - d_M, alphaT_deg);
dCLT_dM = (CL_TpM - CL_TmM) / (2 * d_M);  % per Mach
dCL_dM = (sinA / sinT)*dCLT_dM;

%% === Finite difference for dCm/dM ===
dM = 1e-5;
Cm_TpM = CmT_interp(M0 + d_M, alphaT_deg);
Cm_TmM = CmT_interp(M0 - d_M, alphaT_deg);
dCmT_dM = (Cm_TpM - Cm_TmM) / (2 * d_M);  % per Mach
% Chain rule
dCm_dM = (sinA * cosB / sinT) * dCmT_dM;
% Small Mach step

% Recompute force coefficients at α ± Δα
CD_p = CD_interp(M0 + dM, alphaT_deg);
CL_p = (sinA / sinT) * CLT_interp(M0 + dM, alphaT_deg);
CS_p = (cosA * sinB / sinT) * CLT_interp(M0 + dM, alphaT_deg); % If C_S is also a function of α, interpolate it here

% Recompute force coefficients at α ± Δα
CD_m = CD_interp(M0 - dM, alphaT_deg);
CL_m = (sinA / sinT) * CLT_interp(M0 - dM, alphaT_deg);
CS_m = (cosA * sinB / sinT) * CLT_interp(M0 - dM, alphaT_deg);

% Compute dynamic pressure (assume rho and V0 fixed for this differential)
q_dyn = 0.5 * rho0 * V0^2;

% Forces in body frame at M + dM
[Fx_p, Fy_p, Fz_p] = aeroFrame2apexFrameForces(CD_p * q_dyn * Sref, CL_p * q_dyn * Sref, CS_p * q_dyn * Sref, alpha_rad, beta_rad);
F_p = [Fx_p; Fy_p; Fz_p];
M_p = cross(r_cg, F_p);
Cm_offset_px = Mp(1) / (q_dyn * Sref * cref);
Cm_offset_py = Mp(2) / (q_dyn * Sref * cref);
Cm_offset_pz = Mp(3) / (q_dyn * Sref * cref);

% Forces in body frame at M - dM
[Fx_m, Fy_m, Fz_m] = aeroFrame2apexFrameForces(CD_m * q_dyn * Sref, CL_m * q_dyn * Sref, CS_m * q_dyn * Sref, alpha_rad, beta_rad);
F_m = [Fx_m; Fy_m; Fz_m];
M_m = cross(r_cg, F_m);
Cm_offset_mx = Mm(1) / (q_dyn * Sref * cref);
Cm_offset_my = Mm(2) / (q_dyn * Sref * cref);
Cm_offset_mz = Mm(3) / (q_dyn * Sref * cref);

% Central difference for derivative w.r.t. Mach
dCm_offset_dMx = (Cm_offset_px - Cm_offset_mx) / (2 * dalpha);
dCm_offset_dMy = (Cm_offset_py - Cm_offset_my) / (2 * dalpha);
dCm_offset_dMz = (Cm_offset_pz - Cm_offset_mz) / (2 * dalpha);

% Add CG contribution to total dCm/dM
dCm_dM = dCm_dM + dCm_offset_dMx + dCm_offset_dMy + dCm_offset_dMz;

end