function [epsT0] = Rayleigh(L, x, T01, T02, M, gam)
dT0_dx = (T02-T01)/L;
T0 = dT0_dx*x + T01;

dT0_dx_T0 = dT0_dx/T0;

p1 = (1+gam*M);
p2 = (1+((gam-1)/2*M));
TOIC_T0 = (p1*p2)/(1-M);

epsT0 = TOIC_T0*dT0_dx_T0;

end