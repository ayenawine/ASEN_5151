function [epsA] = area(x, alpha, r1, gam, M)
% conical diffuser
    A = pi*(r1 - x*tand(alpha))^2;
    dA_dx = -2*pi*(r1-x*tand(alpha))*tand(alpha);

    dA_dx_A = dA_dx/A;

    TOIC_A = (-2*(1+((gam-1)/2)*M))/(1-M);

    epsA = TOIC_A*dA_dx_A;
end