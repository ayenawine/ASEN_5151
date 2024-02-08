function [M] = RK4_dMsq(M, xspan, dx, gam, alpha, r1, f, L, T01, T02)


%dMsq = @(x, M, dx, gam, alpha, r1, f, Qdot, mdot, Cp, L, T1) M*((-2*(1+((gam-1)/2)*M))/(1-M)*area(x, alpha, r1)); % + (gam*M(i)*(1+(gam-1)*M(i)))/(1-M(i))*Fanno(xspan(i), alpha, r1, f) + ((1+gam*M(i))*(1+(gam-1)/2)*M(i))/(1-M(i))*Rayleigh(Qdot, mdot, M(i), Cp, L, xspan(i), T1, gam);
%dMsq = @(x, M, dx, gam, alpha, r1, f, Qdot, mdot, Cp, L, T1)  M*((gam*M*(1+((gam-1)/2)*M)/(1-M))*Fanno(x, alpha, r1, f)); % + ((1+gam*M(i))*(1+(gam-1)/2)*M(i))/(1-M(i))*Rayleigh(Qdot, mdot, M(i), Cp, L, xspan(i), T1, gam);
%dMsq = @(x, M, dx, gam, alpha, r1, f, Qdot, mdot, Cp, L, T01, T02) M*(((1+gam*M)*(1+((gam-1)/2)*M))/(1-M))*Rayleigh(Qdot, mdot, Cp, L, x, T01, T02);
%
dMsq = @(x, M, dx, gam, alpha, r1, f, L, T01, T02) M*Rayleigh(L, x, T01, T02, M, gam);
%dMsq = @(x, M, dx, gam, alpha, r1, f, L, T01, T02) M*area(x, alpha, r1, gam, M);
%dMsq = @(x, M, dx, gam, alpha, r1, f, L, T01, T02) M*Fanno(x, alpha, r1, f, gam, M);
%M*(((1+gam*M)*(1+(gam-1)/2)*M)/(1-M)*Rayleigh(Qdot, mdot, M, Cp, L, x, T01, gam));

for i = 1:length(xspan)-1

    m1 = dMsq(xspan(i), M(i), dx, gam, alpha, r1, f, L, T01, T02);
    m2 = dMsq((xspan(i) + dx/2), (M(i) + dx/2*m1), dx, gam, alpha, r1, f, L, T01, T02);
    m3 = dMsq((xspan(i) + dx/2), (M(i) + dx/2*m2), dx, gam, alpha, r1, f, L, T01, T02);
    m4 = dMsq((xspan(i) + dx), (M(i) + dx*m3), dx, gam, alpha, r1, f, L, T01, T02);

    %m1 = (-2*(1+(gam-1)*M(i)))/(1-M(i))*area(xspan(i), alpha, r1) % + (gam*M(i)*(1+(gam-1)*M(i)))/(1-M(i))*Fanno(xspan(i), alpha, r1, f) + ((1+gam*M(i))*(1+(gam-1)/2)*M(i))/(1-M(i))*Rayleigh(Qdot, mdot, M(i), Cp, L, xspan(i), T1, gam);

    %m2 = (-2*(1+(gam-1)*(M(i) + dx/2*m1)))/(1-(M(i) + dx/2*m1))*area((xspan(i) + dx/2), alpha, r1) % + (gam*((M(i)+dx/2*m1)*(1+(gam-1)*(M(i)+dx/2*m1))))/(1-(M(i)+dx/2*m1))*Fanno((xspan(i)+dx/2), alpha, r1, f) + ((1+gam*(M(i)+ dx/2*m1))*(1+(gam-1)/2)*(M(i)+ dx/2*m1))/(1-(M(i)+ dx/2*m1))*Rayleigh(Qdot, mdot, (M(i) + dx/2*m1), Cp, L, (xspan(i)+dx/2), T1, gam);
    %m3 = (-2*(1+(gam-1)*(M(i) + dx/2*m2)))/(1-(M(i) + dx/2*m2))*area((xspan(i) + dx/2), alpha, r1) % + (gam*((M(i)+dx/2*m2)*(1+(gam-1)*(M(i)+dx/2*m2))))/(1-(M(i)+dx/2*m2))*Fanno((xspan(i)+dx/2), alpha, r1, f) + ((1+gam*(M(i)+ dx/2*m2))*(1+(gam-1)/2)*(M(i)+ dx/2*m2))/(1-(M(i)+ dx/2*m2))*Rayleigh(Qdot, mdot, (M(i) + dx/2*m2), Cp, L, (xspan(i)+dx/2), T1, gam);
    %m4 = (-2*(1+(gam-1)*(M(i) + dx/2*m3)))/(1-(M(i) + dx/2*m3))*area((xspan(i) + dx), alpha, r1) % + (gam*((M(i)+dx/2*m3)*(1+(gam-1)*(M(i)+dx/2*m3))))/(1-(M(i)+dx/2*m3))*Fanno((xspan(i)+dx), alpha, r1, f) + ((1+gam*(M(i)+ dx/2*m3))*(1+(gam-1)/2)*(M(i)+ dx/2*m3))/(1-(M(i)+ dx/2*m3))*Rayleigh(Qdot, mdot, (M(i) + dx/2*m3), Cp, L, (xspan(i)+dx), T1, gam);
    

    M(i+1) = M(i) + dx/6*(m1 + 2*m2 + 2*m3 + m4);

end


end