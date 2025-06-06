


% Inputs:
% Q            : 4x1 quaternion [q1; q2; q3; q4] (scalar last)
% V_gamma_chi  : [V; gamma; chi]
% R_tau_lambda : [r; tau; lambda] -> tau: latitude, lambda: longitude

% Extract quaternion elements
q1 = Q(1); q2 = Q(2); q3 = Q(3); q4 = Q(4);

% Direction cosine matrix from body to inertial frame
C_IB = [ ...
    q1^2 + q4^2 - q2^2 - q3^2,   2*(q1*q2 + q3*q4),       2*(q1*q3 - q2*q4);
    2*(q1*q2 - q3*q4),           q2^2 + q4^2 - q1^2 - q3^2, 2*(q2*q3 + q1*q4);
    2*(q1*q3 + q2*q4),           2*(q2*q3 - q1*q4),       q3^2 + q4^2 - q1^2 - q2^2 ...
];

% Extract trajectory angles
gamma = V_gamma_chi(2);  % flight path angle
chi   = V_gamma_chi(3);  % heading angle

% Extract latitude and longitude
tau = R_tau_delta(2);   % latitude
delta = R_tau_delta(3);% longitude

% C_TV = C2(gamma) * C3(chi)
cgamma = cos(gamma);
sgamma = sin(gamma);
cchi = cos(chi);
schi = sin(chi);
C_TV = [cgamma*cchi,  schi*cgamma,  -sgamma;
       -schi,     cchi,      0;
        cchi*sgamma,  schi*sgamma,   cgamma];

% C_VI = C2(-delta - pi/2) * C3(tau)
cdelta = cos(delta);
sdelta = sin(delta);
ct = cos(tau);
st = sin(tau);

C_VI = [ ...
    -ct*sdelta,  -st*sdelta,  cdelta;
    -st,      ct,     0;
    -ct*cdelta,  -st*cdelta, -sdelta
];

% Compose full transformation: C_TB
C_TB = C_TV * C_VI * C_IB;

% Extract angles
alpha  = atan2(C_TB(1,3), C_TB(1,1));   % [rad]
beta   = asin(C_TB(1,2));               % [rad]
sigma  = atan2(-C_TB(3,2), C_TB(2,2));  % [rad]

% Output
alpha_beta_sigma = [alpha; beta; sigma];
