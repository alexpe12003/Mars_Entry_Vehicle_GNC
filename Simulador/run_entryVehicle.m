% =====================================================================
% Initialization Script for Entry Vehicle Simulation (Earth Case)
% =====================================================================
%
% This script sets up all necessary constants, initial conditions, and
% vehicle/aerodynamic properties for simulating the atmospheric entry of a
% capsule-like vehicle (Apollo-class) into Earth's atmosphere.
%
% It includes:
%   - Definition of planetary and atmospheric parameters (Earth)
%   - Loading of aerodynamic coefficient datasets (drag, lift, moments)
%   - Specification of simulation time settings
%   - Initialization of the vehicle's translational and rotational state
%   - Computation of initial attitude quaternion and ECI coordinates
%   - Definition of physical properties such as mass, inertia, and reference geometry
%   - Configuration of inertial measurement unit (IMU) sensor biases and errors
%
% The initialized variables are later used in the main Simulink model and
% associated control/guidance scripts.
%
% Authors:
%   Alexandre Alves Pereira       Nº112288
%   Eduardo De Almeida Helena     Nº102793
%   Guilherme Alexandre Pires Martins Nº112118
%   Rafael Coimbra Azeiteiro      Nº102478
%   Rodrigo Marques Pereira       Nº112117
%
% Professor:
%   Alain de Souza
%
% Course Project — Guidance, Navigation and Control, IST Lisbon University
% =====================================================================

clc;
clear all;

%% Planet Constants - Earth

% === Earth Radius ===
Mars_radius = 6371e3;  % Radius of Earth in meters [m]

% === Earth Gravitational Parameter ===
Earth_mu = 3.98319648e14;  % Gravitational parameter for Earth [m^3/s^2]

% === Earth Atmospheric Parameters ===
H_s = 7050;  % Scale height of Earth' atmosphere [m]
rho_0 = 1.225;  % Reference density of Earth' atmosphere at sea level [kg/m^3]

% === Specific Heat Ratio ===
specific_heat = 1.29;  % Specific heat ratio for Earth' atmosphere (Cp/Cv), dimensionless

% === Specific Gas Constant for Earth ===
R_star = 8314.32;  % Universal gas constant [J/(kmol*K)]
M_CO2 = 44.01;  % Molecular mass of CO2 [kg/kmol]
R_s = 287;  % Specific gas constant for Earth' atmosphere [J/(kg*K)]


%% Load Aerodynamic Data
load("ApolloDragCoeff.mat");
load("ApolloLiftCoeff.mat");
load("ApolloPitchMomentCoeff.mat");
load("ApolloPitchDampingCoeff.mat");


%% Simulation Settings
startTime = 0.0;  % Start time of simulation [s]
stopTime = 1200;  % Stop time of simulation [s]
stepTime = 0.01;  % Time step of simulation [s]

%% Initial Conditions

% === Translation State ===
V_0 = 11e3;  % Initial velocity (ground speed) [m/s]
gamma_0 = deg2rad(-9.536);  % Initial flight-path angle [rad]
chi_0 = deg2rad(90);  % Initial heading angle [rad]

R_0 = 220e3 + Mars_radius;  % Initial position radius [m] (distance from Mars center)
tau_0 = deg2rad(0);  % Initial longitude [rad]
lambda_0 = deg2rad(0);  % Initial latitude [rad]

% === Attitude State ===
p_0 = deg2rad(0);  % Initial roll rate [rad/s]
q_0 = deg2rad(0);  % Initial pitch rate [rad/s]
r_0 = deg2rad(0);  % Initial yaw rate [rad/s]

alpha_0 = deg2rad(-23.82);  % Initial angle of attack [rad]
beta_0 = deg2rad(0.0);  % Initial angle of sideslip [rad]
sigma_0 = deg2rad(100.);  % Initial bank angle [rad]

% === Initial Mass ===
m_0 = 4976;  % Initial mass of the vehicle [kg]

% === Initial Moments of Inertia ===
I_xx = 5617.61;  % Moment of inertia about the x-axis [kg/m^2]
I_yy = 4454.62;  % Moment of inertia about the y-axis [kg/m^2]
I_zz = 4454.80;  % Moment of inertia about the z-axis [kg/m^2]

% === Initial Products of Inertia ===
I_xy = 0;  % Product of inertia in the xy-plane [kg/m^2]
I_yz = 0;  % Product of inertia in the yz-plane [kg/m^2]
I_xz = 0;  % Product of inertia in the xz-plane [kg/m^2]

% Get initial quaternion values and initial ECI coordinates
Q0 = get_Q0(tau_0, lambda_0, chi_0, gamma_0, alpha_0, beta_0, sigma_0);  % Compute initial quaternion (orientation)
[r0_eci, v0_eci] = getInitialECI(R_0, tau_0, lambda_0, V_0, gamma_0, chi_0);  % Initial position and velocity in ECI (Earth-Centered Inertial)


%% Vehicle Properties

% === Reference Area ===
S_ref = 12;  % Reference area of the vehicle [m^2]

% === Reference Length ===
d_ref = 3.9;  % Reference diameter of the vehicle [m]

% === Center of Mass Location ===
r_cm = [-0.137; 0; 1.8];  % Center of mass location relative to the apex [m]


%% IMU Sensor Properties

% === Accelerometer Properties ===
acce_bias = [0; 0; 0];  % Accelerometer bias [m/s^2]
acce_S = [0, 0, 0;  % Accelerometer sensitivity matrix
          0, 0, 0;
          0, 0, 0];

% === Rate Gyroscope Properties ===
rateGyro_bias = [0; 0; 0];  % Gyroscope bias [rad/s]
rateGyro_S = [0, 0, 0;  % Gyroscope sensitivity matrix
              0, 0, 0;
              0, 0, 0];
