function dM2_dx = dp_combined_M(x,M,constants)

    dM2_dx_area = dp_area_M(x,M,constants);
    dM2_dx_fanno = dp_fanno_M(x,M,constants);
    dM2_dx_rayleigh = dp_rayleigh_M(x,M,constants);

    dM2_dx = dM2_dx_area + dM2_dx_fanno + dM2_dx_rayleigh;
end