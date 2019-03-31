function [Fvec] = compute_loads(t, H, x4, x5)

omega   = H.omega;
zeta    = H.zeta;
k       = H.k;
rho     = H.rho;
g       = H.g;

Bpp     = H.Bpp;
Mpp     = H.Mpp;
Tpp     = H.Tpp;
Lpp     = H.Lpp;
B_0     = H.B_0;
a       = H.a;
b       = H.b;
Mb      = H.Mb;
Ml      = H.Ml;

Fh      = (-rho*g*zeta*cos(1/2*k*Bpp+omega*t)*(exp(k*Tpp)*Tpp*k-exp(k*Tpp)+1)/k^2/Tpp+rho*g*zeta*cos(1/2*k*Bpp-omega*t)*(exp(k*Tpp)*Tpp*k-exp(k*Tpp)+1)/k^2/Tpp)*Lpp;
Fv      = (-exp(-k*Tpp)*g*rho*zeta*(sin(1/2*k*Bpp+omega*t)*k*Bpp-2*cos(omega*t)+2*cos(1/2*k*Bpp+omega*t))/Bpp/k^2+exp(-k*Tpp)*g*rho*zeta*(sin(1/2*k*Bpp-omega*t)*k*Bpp+2*cos(1/2*k*Bpp-omega*t)-2*cos(omega*t))/Bpp/k^2)*Lpp;
Mwaves  = Lpp*(-1/2*exp(-k*Tpp)*g*rho*zeta*(sin(1/2*k*Bpp+omega*t)*k*Bpp-2*cos(omega*t)+2*cos(1/2*k*Bpp+omega*t))/k^2-1/2*exp(-k*Tpp)*g*rho*zeta*(sin(1/2*k*Bpp-omega*t)*k*Bpp+2*cos(1/2*k*Bpp-omega*t)-2*cos(omega*t))/k^2+rho*g*zeta*cos(1/2*k*Bpp+omega*t)*(exp(k*Tpp)*Tpp*k-exp(k*Tpp)+1)/k^2-rho*g*zeta*cos(1/2*k*Bpp-omega*t)*(exp(k*Tpp)*Tpp*k-exp(k*Tpp)+1)/k^2);
Mmasses = Mb*g*x5+Ml*g*x4;
Fmasses = H.Fmasses; % correct?

Mtheta2 = 0;
Mtheta3 = 0;

Fvec    = [Fh; Fv + Fmasses; Mwaves + Mmasses; Mtheta2; Mtheta3;];

end