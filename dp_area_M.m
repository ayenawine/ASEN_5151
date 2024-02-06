function dM2_dx = dp_area_M(x,M_in,constants)
    p1 = (1+((constants.gamma-1)/2)*M_in^2);
    dAdx_A = ( 2 * tan(constants.alpha) ) / ( constants.r_1 + x*tan(constants.alpha) );
    
    dM2_dx = ( ( -2*p1 ) / (1-M_in^2) )  *  dAdx_A * M_in^2;
end