% =====================================================================
% Entry Vehicle Simulation Setup Script
% =====================================================================
%
% This script serves as the main driver for initializing and running all 
% major components involved in the entry vehicle simulation. It is intended
% to be executed before performing further analysis (e.g., state-space modeling,
% guidance computation, Monte Carlo runs).
%
% The sequence of operations includes:
%   1. Opening the Simulink model for the entry vehicle ('entryVehicle.slx')
%   2. Running the base simulation script to initialize nominal trajectory
%   3. Executing the state-space modeling script
%   4. (Optional) Running the navigation script
%   5. (Optional) Running the guidance script to generate nominal trajectory data
%
% Notes:
%   - The guidance script, when enabled, retrieves and saves the nominal 
%     trajectory described in the report to the 'data_sets_guidance' directory.
%   - Ensure all relevant `.m` and `.slx` files are in the working directory.
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
run('guidance.m');        % Execute the guidance.m script
% Upon running the guidance script, the program retrieves the information
% about the nominal trajectory described in the report, currently stored
% in the 'data_sets_guidance' directory.
