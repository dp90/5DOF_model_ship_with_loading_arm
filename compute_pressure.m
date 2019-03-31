function p = compute_pressure(x, z, t, H)

p = H.rho*H.g*H.zeta*exp(-H.k*z)*cos(H.k*x - H.omega*t);

end

