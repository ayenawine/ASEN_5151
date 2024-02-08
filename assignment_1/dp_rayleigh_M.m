function dM2_dx = dp_rayleigh_M(x,M,constants)

    p1 = (1 + constants.gamma*M);
    p2 = ( 1 + ((constants.gamma-1)/2)*M );

    T_0 = constants.T_01 + constants.dT_0_dx * x;

    dM2_dx = ( (p1*p2) / (1-M) )  * M * (constants.dT_0_dx/T_0);
end