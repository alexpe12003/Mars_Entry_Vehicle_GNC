function [Q0] = get_Q0(tau_0, lambda_0, chi_0, gamma_0, alpha_0, beta_0, sigma_0)
   
    % Precompute sines and cosines
    ct = cos(tau_0);
    st = sin(tau_0);
    cd = cos(lambda_0);
    sd = sin(lambda_0);
    
    cx = cos(chi_0);
    sx = sin(chi_0);
    cg = cos(gamma_0);
    sg = sin(gamma_0);
    
    cs = cos(sigma_0);
    ss = sin(sigma_0);
    cb = cos(beta_0);
    sb = sin(beta_0);
    ca = cos(alpha_0);
    sa = sin(alpha_0);
    
    % Matrix C_I,V
    C_IV = [ -ct*sd,      -st,   -ct*cd;
             -st*sd,       ct,   -st*cd;
                 cd,        0,       -sd ];
    
    % Matrix C_V,T
    C_VT = [ cx*cg,     -sx,   cx*sg;
             sx*cg,      cx,   sx*sg;
               -sg,       0,      cg ];
    
    % Matrix C_T,B
    C_TB = [ ca*cb,         sb,               sa*cb;
            -ca*sb*cs - sa*ss, cb*cs, ca*ss - sa*sb*cs;
             ca*sb*ss - sa*cs, -cb*ss, ca*cs + sa*sb*ss ];

    C_IB = C_IV * C_VT * C_TB;
    
    % Compute Q4 first
    Q4 = 0.5 * sqrt(1 + C_IB(1,1) + C_IB(2,2) + C_IB(3,3));
    
    % Compute Q1, Q2, Q3 based on Q4
    Q1 = (C_IB(2,3) - C_IB(3,2)) / (4 * Q4);
    Q2 = (C_IB(3,1) - C_IB(1,3)) / (4 * Q4);
    Q3 = (C_IB(1,2) - C_IB(2,1)) / (4 * Q4);
    
    % Quaternion as a column vector (or row vector if preferred)
    Q0 = [Q1; Q2; Q3; Q4];
    Q0 = Q0/norm(Q0);
end
