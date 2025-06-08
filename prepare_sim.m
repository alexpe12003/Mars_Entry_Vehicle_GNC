clc;clear
% === Open Entry Vehicle Simulink Model ===
open_system('entryVehicle.slx');  % Open the Simulink model

% === Run Entry Vehicle Simulation ===
run('run_entryVehicle.m');  % Run the simulation

% === Run State Space Model ===
run('stateSpace.m');      % Execute the stateSpace.m script

% === Run Navigation Script ===
%run('navigation.m');        % Execute the navigation.m script

% === Run Guidance Script ===
%run('guidance.m');        % Execute the guidance.m script
