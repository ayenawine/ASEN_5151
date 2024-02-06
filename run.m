% cleanup
clear
clc
close all

% constant inputs
gamma = 1.4;
alpha = deg2rad(1); % radians
R = 287; % J / kg * K
f_f = 0.025; % m
Q_dot = 500*1000; % J/s
L = 8.0; % m
dt = 0.25; % m

% inlet conditions
r_1 = 0.5; % m
M_1 = 0.2;

xs = 0:dt:L;
Ms = zeros([1 length(xs)]);
Ms(1) = M_1;

for i = 1:length(xs)-1
    x = xs(i);
    M = Ms(i);

    xn = x;
    yn = M;
    m1 = dp_area(xn,yn,gamma,alpha,r_1);

    xn = x+dt/2;
    yn = M+(dt/2)*m1;
    m2 = dp_area_M(xn,yn,gamma,alpha,r_1);

    xn = x+dt/2;
    yn = M+(dt/2)*m2;
    m3 = dp_area_M(xn,yn,gamma,alpha,r_1);

    xn = x+dt;
    yn = M+dt*m3;
    m4 = dp_area_M(xn,yn,gamma,alpha,r_1);

    Ms(i+1) = M + (dt*(m1+2*m2+2*m3+m4)) / 6;
end

figure
plot(xs,Ms)
xlabel("Duct Length [m]")
ylabel("Mach")
grid