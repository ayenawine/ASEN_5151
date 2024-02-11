%% Notes
% CONVERGING DUCT
% 
% Run script for testing driving potential functions against table for problem 2
%

%%

% cleanup
clear
clc
close all

% gather all inputs to a struct so we can share for all potential functions
constants = {};

% key input variables
constants.dx = 0.25; % m
constants.M_1 = 0.2;

% constant inputs
constants.gamma = 1.4;
constants.alpha = deg2rad(1); % radians
constants.R = 287; % J / kg * K
constants.f_f = 0.025; % m
constants.Q_dot = 500*1000; % J / s
constants.L = 8.0; % m
constants.c_p = (constants.gamma*constants.R) / (constants.gamma - 1); % section 1, page 8, cp constant for certain temperature range (ones we care about)

% simulation flags
constants.converging = true;
constants.combined_potential = false;

% inlet conditions
constants.r_1 = 0.5; % m, radius
constants.A_1 = pi*constants.r_1^2; % m^2, area
constants.T_1 = 300; % K, temperature
constants.P_1 = 100000; % Pa, kg/m*s2, pressure
constants.rho_1 = constants.P_1 / (constants.R * constants.T_1); % kg/m3, density, section 1, page 6
constants.V_1 = constants.M_1 * sqrt(constants.gamma*constants.R*constants.T_1); % m/s, Section 1, page 10, or https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/mchkdrv.html
constants.m_dot1 = constants.rho_1 * constants.V_1 * constants.A_1; % kg / s

% calculate T_0 at position 1 and 2 for later
constants.T_01 = (1 + ((constants.gamma-1)/2)*constants.M_1^2) * constants.T_1; % K, section 1, page 13
constants.T_02 = constants.T_01 + constants.Q_dot / (constants.m_dot1 * constants.c_p); % K, energy equation, section 2, page 5

% from problem statement, dT_0/dx is assumed constant from heating so we can do this
constants.dT_0_dx = (constants.T_02 - constants.T_01) / constants.L; % K / m

xs = 0:constants.dx:constants.L;
M2_area = zeros([1 length(xs)]);
M2_fanno = zeros([1 length(xs)]);
M2_rayleigh = zeros([1 length(xs)]);
M2_area(1) = constants.M_1^2;
M2_fanno(1) = constants.M_1^2;
M2_rayleigh(1) = constants.M_1^2;

for i = 1:length(xs)-1
    x = xs(i);

    % NOTE TO SELF:
    % ODE is M2 with respect to dx, so M2 has to be input into the runge
    % kutta, not M

    % calculate Area driven potential using runge kutta 4th order
    M2_area(i+1) = myRungeKutta4(@dp_area_M,constants,M2_area(i),x,constants.dx);

    % calculate Friction driven potential using runge kutta 4th order
    M2_fanno(i+1) = myRungeKutta4(@dp_fanno_M,constants,M2_fanno(i),x,constants.dx);

    % calculate Heat driven potential using runge kutta 4th order
    M2_rayleigh(i+1) = myRungeKutta4(@dp_rayleigh_M,constants,M2_rayleigh(i),x,constants.dx);
end

figure
hold on 
plot(xs,sqrt(M2_area),"DisplayName",sprintf("Area Driven, M_L = %.5f",sqrt(M2_area(end))))
plot(xs,sqrt(M2_fanno),"DisplayName",sprintf("Friction Driven, M_L = %.5f",sqrt(M2_fanno(end))))
plot(xs,sqrt(M2_rayleigh),"DisplayName",sprintf("Heat Driven, M_L = %.5f",sqrt(M2_rayleigh(end))))
legend("Location","best")
title("Mach, individual potentials")
xlabel("Duct Length [m]")
ylabel("Mach")
grid

%% calculate T/T1 for each driving potential

% calculate stagnation temperature across the duct using dT0_dx (rayleigh flow only)
T_0x = constants.T_01 + constants.dT_0_dx * xs;

% for area change or friction only, assume no heat addition
T_02 = constants.T_01;

% calculate temperature distribution, area driven
num = 1 + ((constants.gamma-1)/2)*constants.M_1^2; % scalar
den = 1 + ((constants.gamma-1)/2).*M2_area; % array
Tx_T1_area = (T_02/constants.T_01).*(num./den); % isentropic relation, section 2, page 32
Tx_area = Tx_T1_area.*constants.T_1;

% calculate temperature distribution, friction driven
num = 1 + ((constants.gamma-1)/2)*constants.M_1^2; % scalar
den = 1 + ((constants.gamma-1)/2).*M2_fanno; % array
Tx_T1_fan = (T_02/constants.T_01).*(num./den); % isentropic relation, section 2, page 32
Tx_fan = Tx_T1_fan.*constants.T_1;

% calculate temperature distribution, heat driven
num = 1 + ((constants.gamma-1)/2)*constants.M_1^2; % scalar
den = 1 + ((constants.gamma-1)/2).*M2_rayleigh; % array
Tx_T1_ray = (T_0x/constants.T_01).*(num./den); % isentropic relation, section 2, page 32
Tx_ray = Tx_T1_ray.*constants.T_1;

figure
hold on 
plot(xs,Tx_T1_area,"DisplayName",sprintf("Area Driven, Tx/T1_L = %.5f",Tx_T1_area(end)))
plot(xs,Tx_T1_fan, "DisplayName",sprintf("Friction Driven, Tx/T1_L = %.5f",Tx_T1_fan(end)))
plot(xs,Tx_T1_ray, "DisplayName",sprintf("Heat Driven, Tx/T1_L = %.5f",Tx_T1_ray(end)))
legend("Location","best")
title("T/T1, individual potentials")
xlabel("Duct Length [m]")
ylabel("T/T1")
grid

%% calculate P/P1 for area driving potential

% calculate velocity distribution
Mx_area = sqrt(M2_area);
Vx_V1_area = (Mx_area/constants.M_1).*sqrt(Tx_T1_area); % definition of Mach number, section 2, page 33
%Vx = Vx_V1.*constants.V_1;

% calculate density distribution
Ax = pi * (constants.r_1-xs*tan(constants.alpha)).^2; % for converging duct
A1_Ax = constants.A_1 ./ Ax;
rhox_rho1_area = A1_Ax ./ (Vx_V1_area); % continuity equation, section 2, page 33

% calculate pressure distribution
Px_P1_area = rhox_rho1_area ./ Tx_T1_area; % equation of state, section 2, page 33
%Px_area = Px_P1_area * constants.P_1;

%% calculate P/P1 for friction driving potential

% calculate velocity distribution
Mx_fan = sqrt(M2_fanno);
Vx_V1_fan = (Mx_fan/constants.M_1).*sqrt(Tx_T1_fan); % definition of Mach number, section 2, page 33
%Vx = Vx_V1.*constants.V_1;

% calculate density distribution
% assume area change is zero
rhox_rho1_fan = 1 ./ (Vx_V1_fan); % continuity equation, section 2, page 33

% calculate pressure distribution
Px_P1_fan = rhox_rho1_fan ./ Tx_T1_fan; % equation of state, section 2, page 33
%Px_fan = Px_P1_fan * constants.P_1;

%% calculate P/P1 for heat driving potential

% calculate velocity distribution
Mx_ray = sqrt(M2_rayleigh);
Vx_V1_ray = (Mx_ray/constants.M_1).*sqrt(Tx_T1_ray); % definition of Mach number, section 2, page 33
%Vx = Vx_V1.*constants.V_1;

% calculate density distribution
% assume area change is zero
rhox_rho1_ray = 1 ./ (Vx_V1_ray); % continuity equation, section 2, page 33

% calculate pressure distribution
Px_P1_ray = rhox_rho1_ray ./ Tx_T1_ray; % equation of state, section 2, page 33
%Px_fan = Px_P1_ray * constants.P_1;

% plot
figure
hold on 
plot(xs,Px_P1_area,"DisplayName",sprintf("Area Driven, Px/P1_L = %.5f",Px_P1_area(end)))
plot(xs,Px_P1_fan,"DisplayName",sprintf("Friction Driven, Px/P1_L = %.5f",Px_P1_fan(end)))
plot(xs,Px_P1_ray,"DisplayName",sprintf("Heat Driven, Px/P1_L = %.5f",Px_P1_ray(end)))
legend("Location","best")
title("P/P1, individual potentials")
xlabel("Duct Length [m]")
ylabel("P/P1")
grid