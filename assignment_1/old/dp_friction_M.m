function dM2_dx = dp_friction_M(x,M,constants)
    p1 = (1+((constants.gamma-1)/2)*M);

    if constants.combined_potential == true
        if constants.converging == true
            D_h = 2*(constants.r_1-x*tan(constants.alpha));
        else
            D_h = 2*(constants.r_1+x*tan(constants.alpha));
        end
    else
        % for constant area assumption
        D_h = constants.r_1 * 2;
    end

    dM2_dx = ( ( constants.gamma*p1*M ) / (1-M) ) * M * (constants.f_f/D_h);
end