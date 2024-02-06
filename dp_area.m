function M_out = dp_area_M(x,M_in,gamma,alpha,r_1)
    M2_in = M_in^2;
    dAdx_A = ( 2 * tan(alpha) ) / ( r_1 + x*tan(alpha) );
    
    M_out = ( ( -2*(1+((gamma-1)/2)*M2_in) ) / (1-M2_in) ) * dAdx_A * M2_in;
end