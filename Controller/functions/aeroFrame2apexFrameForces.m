function [X_A, Y_A, Z_A] = aeroFrame2apexFrameForces(D, L, S, alpha, beta)

c_alpha = cos(alpha);
s_alpha = sin(alpha);
c_beta = cos(beta);
s_beta = sin(beta);

X_A = -D*c_alpha*c_beta + S*c_alpha*s_beta + L*s_alpha;
Y_A = -D*s_beta - S*c_beta;
Z_A = -D*s_alpha*c_beta + S*s_alpha*s_beta - L*c_alpha;
end
