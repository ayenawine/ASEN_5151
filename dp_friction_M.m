function dM2_dx = dp_friction_M(x,M_in,constants)
    p1 = (1+((constants.gamma-1)/2)*M_in^2);
    D_h = 2*(constants.r_1+x*tan(constants.alpha));

    dM2_dx = ( ( constants.gamma*p1*M_in^2 ) / (1-M_in^2) )  *  M_in^2 * (constants.f_f/D_h);
end