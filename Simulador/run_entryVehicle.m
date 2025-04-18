
% Load Aerodynamic Data
load("ApolloDragCoeff.mat");
load("ApolloLiftCoeff.mat");
load("ApolloPitchMomentCoeff.mat");
load("ApolloPitchDampingCoeff.mat");

% Constants
Mars_radius = 3389.5e3; % [m]
Mars_mu = 4.27985580e13; % [m^3/s^2]

% Mars atmosphere
H_s = 8.8e3; % [m] - Scale height
rho_0 = 0.0168; % [kg/m^3] - Reference density

specific_heat = 1.29; % [-] - Specific heat ratio

R_star = 8314.32; % [J/(kmol*K)] - Universal gas constant
M_CO2 = 44.01; % [kg/kmol] - CO2 molecular mass
R_s = R_star/M_CO2; % [J/(kg*K)] - Specific gas constant


% Simulation
startTime = 0.0; % [s]
stopTime = 200.0; % [s]

stepTime = 0.0001; % [s]

%% Initial Conditions

% Translation State
V_0 = 11e3; % [m/s] - Ground speed
gamma_0 = deg2rad(-9.536); % [rad] - Flight-path angle
chi_0 = deg2rad(90); % [rad] - Heading

R_0 = 70e3 + Mars_radius; % [m] - Position radius
tau_0 = deg2rad(0); % [rad] - Longitude
lambda_0 = deg2rad(0); % [rad] - Latitude

% Attitude State
p_0 = deg2rad(0); % [rad/s] - Roll rate
q_0 = deg2rad(0); % [rad/s] - Pitch rate
r_0 = deg2rad(0); % [rad/s] - Yaw rate

alpha_0 = deg2rad(-24.5); % [rad] - Angle of attack
beta_0 = deg2rad(0.05); % [rad] - Angle of sideslip
sigma_0 = deg2rad(180); % [rad] - Bank angle

% Initial mass
m_0 = 4976; % [kg]

% Initial moments of inertia
I_xx = 5617.61; % [kg/m^2]
I_yy = 4454.62; % [kg/m^2]
I_zz = 4454.80; % [kg/m^2]

% Initial products of inertia
I_xy = 0; % [kg/m^2]
I_yz = 0; % [kg/m^2]
I_xz = 1752.00; % [kg/m^2]


%% Vehicle Properties

% Reference Area
S_ref = 12; % [m^2]

% Reference Length
d_ref = 3.9; % [m^2] - Vehicle diameter
