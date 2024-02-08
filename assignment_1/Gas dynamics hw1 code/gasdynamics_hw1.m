%% 
% ASEN 5151: Gas Dynamics 
% Spring 2024
% Assignment 1
% Due: 2/9, 5 PM
% Cate Leszcz

%% Problem Setup:
gam = 1.4;
R = 287; % J/kgK
r1 = 50e-2; % m
alpha = 0; % deg
f = 0.025;
Qdot = 500e3; % J/s
L = 8; % m

% initial conditions
T1 = 300; % K
P1 = 1e5; % Pa
M1 = [0.2^2] %, 3];
dx = [0.25]; %, 0.5, 1]; % m

% calculating needed values
Cp = (gam*R)/(gam-1);
rho1 = P1/(R*T1);
V1 = sqrt(M1)*sqrt(gam*R*T1);
A1 = pi*r1^2;
mdot = rho1*V1*A1;
T01 = T1*(1 +(gam-1)/2*M1);
T02 = Qdot/(mdot*Cp) + T01;

xspan = 0:dx:L;

% calling Runge-Kutta function
%Msq = RK4(M1(1), xspan, dx(1), gam, r1, alpha)
M = RK4_dMsq(M1, xspan, dx, gam, alpha, r1, f, L, T01, T02);
%Msq = RK4(M1, xspan, dx, gam, r1, alpha, f, mdot, Cp, L, Qdot, T0)


% solving for M
M2 = sqrt(M(end))

% plot(xspan, sqrt(M))

