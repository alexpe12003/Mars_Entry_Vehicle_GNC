% Create Apollo pitch-moment coefficient lookup table data for Simulink
% Based on Table B.3 from NAA 1965

% Mach numbers (first column in the table)
Mach_numbers_Cm = [0.90; 1.10; 1.20; 1.35; 1.65; 2.00; 2.40; 3.00; 4.00; 6.00; 10.00];

% Angles of attack (first row in the table)
alpha_T_Cm = [0, 5, 10, 15, 20, 25, 30, 35, 40];

% Pitch-moment coefficient data (Cm_T)
Cm_T = [...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0944, -0.1690, -0.2479, -0.3098, -0.3546, -0.3888, -0.4188, -0.4402;...
    0.0, -0.0902, -0.1754, -0.2501, -0.3226, -0.3695, -0.3974, -0.4188, -0.4360;...
    0.0, -0.0951, -0.1899, -0.2794, -0.3400, -0.4006, -0.4322, -0.4586, -0.4692;...
    0.0, -0.0951, -0.1846, -0.2820, -0.3505, -0.4111, -0.4428, -0.4718, -0.4797];

% Save the data to a .mat file
save('ApolloPitchMomentCoeff.mat', 'Mach_numbers_Cm', 'alpha_T_Cm', 'Cm_T');

% Display confirmation
disp('Apollo pitch-moment coefficient data saved to ApolloPitchMomentCoeff.mat');
disp('To use in Simulink:');
disp('1. Load the .mat file in your model workspace');
disp('2. Use a 2D Lookup Table block');
disp('3. Set the Row index to "Mach_numbers"');
disp('4. Set the Column index to "alpha_T"');
disp('5. Set the Table data to "Cm_T"');
disp(' ');
disp('Note: Connect Mach number to u1 (top input) and angle of attack to u2 (bottom input)');
disp('Negative values indicate nose-down pitching moments');