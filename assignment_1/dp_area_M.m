function dM2_dx = dp_area_M(x,M,constants)
    % Example: Section 2, page 34
    p1 = ( 1 + ((constants.gamma-1)/2)*M );

    if constants.converging == true
        dAdx_A = -( 2 * tan(constants.alpha) ) / ( constants.r_1 - x*tan(constants.alpha) );
    else
        dAdx_A = ( 2 * tan(constants.alpha) ) / ( constants.r_1 + x*tan(constants.alpha) );
    end
    
    dM2_dx = ( ( -2*p1 ) / (1-M) )  *  dAdx_A * M;
end