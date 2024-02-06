function Mout = myRungeKutta4(dp_ODE,constants,M,x,dx)

    xn = x;
    yn = M;
    m1 = dp_ODE(xn,yn,constants);

    xn = x+dx/2;
    yn = M+(dx/2)*m1;
    m2 = dp_ODE(xn,yn,constants);

    xn = x+dx/2;
    yn = M+(dx/2)*m2;
    m3 = dp_ODE(xn,yn,constants);

    xn = x+dx;
    yn = M+dx*m3;
    m4 = dp_ODE(xn,yn,constants);

    Mout = M + ( dx*( m1 + 2*m2 + 2*m3 + m4 ) ) / 6;
end