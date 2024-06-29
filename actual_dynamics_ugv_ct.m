function dx_dt = actual_dynamics_ugv_ct(t, x, u)
    % continous-time update equation of the nominal dynamics in earth frame.
    % For a unicycle
    dxdt(1) = u(1)*cos(x(3));
    dxdt(2) = u(1)*sin(x(3));
    dxdt(3) = u(2);

end