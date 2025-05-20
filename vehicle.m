classdef vehicle
    properties
        m_0
        I_xx
        I_yy
        I_zz
        I_xy
        I_yz
        I_xz
        S_ref
        d_ref
    end

    methods

        function obj = parafoil(m_0, I_xx, I_yy, I_zz, I_xy, I_yz, I_xz)
            obj.m0 = m_0;
            obj.Ixx = I_xx;
            obj.Iyy = I_yy;
            obj.Izz = I_zz;
            obj.Ixy = I_xy;
            obj.Iyz = I_yz;
            obj.Ixz = I_xz;
            obj.S_ref = S_ref;
            obj.d_ref = d_ref;
        end
    end
end