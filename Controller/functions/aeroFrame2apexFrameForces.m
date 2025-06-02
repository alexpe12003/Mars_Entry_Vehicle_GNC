function [X_A, Y_A, Z_A] = aeroFrame2apexFrameForces(D, L, S, alpha, beta)
% Transform aerodynamic forces (Drag D, Lift L, Side force S)
% from the aerodynamic frame to the apex (body-fixed) frame.
% 
% Inputs:
%   D     - Drag force (along the freestream direction)
%   L     - Lift force (perpendicular to the freestream in the vertical plane)
%   S     - Side force (perpendicular to the freestream in the horizontal plane)
%   alpha - Angle of attack (radians)
%   beta  - Sideslip angle (radians)
%
% Outputs:
%   X_A - Force component along the apex x-axis (body-forward)
%   Y_A - Force component along the apex y-axis (body-right)
%   Z_A - Force component along the apex z-axis (body-down)


% Precompute sine and cosine of angles for efficiency
c_alpha = cos(alpha);
s_alpha = sin(alpha);
c_beta = cos(beta);
s_beta = sin(beta);

% Transform forces from aerodynamic to apex frame using standard rotation:
% First rotate around the z-axis by -beta, then around the y-axis by -alpha.
X_A = -D*c_alpha*c_beta + S*c_alpha*s_beta + L*s_alpha;
Y_A = -D*s_beta - S*c_beta;
Z_A = -D*s_alpha*c_beta + S*s_alpha*s_beta - L*c_alpha;
end
