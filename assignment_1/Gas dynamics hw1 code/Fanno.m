function [epsf] = Fanno(x, alpha, r1, f, gam, M)
    D_h = 2*(r1 - x*tand(alpha))

    f_Dh = f/D_h;

    TOIC_f = (gam*M*(1+((gam-1)/2)*M)/(1-M));

    epsf = TOIC_f*f_Dh;
end