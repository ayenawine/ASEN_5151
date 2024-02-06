function dM2_dx = dp_heat_M(x,M_in,constants)
    p1 = (1 + constants.gamma*M_in^2);
    p2 = ( 1 + ((constants.gamma-1)/2)*M_in^2 );

    T_0 = constants.T_01 + constants.dT_0_dx * x;

    dM2_dx = (( p1*p2 ) / (1-M_in^2))  *  M_in^2 * (constants.dT_0_dx/T_0);
end