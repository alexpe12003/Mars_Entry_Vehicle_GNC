function [r0_eci, v0_eci] = getInitialECI(R_0, tau_0, lambda_0, V_0, gamma_0, chi_0)

% Position in ECI
x = R_0 * cos(lambda_0) * cos(tau_0);
y = R_0 * cos(lambda_0) * sin(tau_0);
z = R_0 * sin(lambda_0);
r0_eci = [x; y; z];

% Velocity in NED (North-East-Down)
vN = V_0 * cos(gamma_0) * cos(chi_0);
vE = V_0 * cos(gamma_0) * sin(chi_0);
vD =  -V_0 * sin(gamma_0);
v_ned = [vN; vE; vD];

% Transformation matrix from NED to ECI
sL = sin(lambda_0);
cL = cos(lambda_0);
sT = sin(tau_0);
cT = cos(tau_0);
T_ECI_NED = [
    -sL*cT, -sT, -cL*cT;
    -sL*sT,  cT, -cL*sT;
     cL,     0, -sL
];

% Velocity in ECI
v0_eci = T_ECI_NED * v_ned;
end
