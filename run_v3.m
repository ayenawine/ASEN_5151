%% Notes
% should have 6 total output plots

%%

% cleanup
clear
clc
close all

% constant inputs
constants = {};
constants.gamma = 1.4;
constants.alpha = deg2rad(1); % radians
constants.R = 287; % J / kg * K
constants.f_f = 0.025; % m
constants.Q_dot = 500*1000; % J / s
constants.L = 8.0; % m
constants.dx = 0.25; % m
constants.c_p = (constants.gamma*constants.R) / (constants.gamma - 1); % section 1, page 8, cp constant for certain temperature range (ones we care about)

% inlet conditions
constants.r_1 = 0.5; % m
constants.A_1 = pi*constants.r_1^2; % m^2
constants.M_1 = 0.2;
constants.T_1 = 300; % K
constants.P_1 = 10000; % Pa, kg/m*s2
constants.rho_1 = constants.P_1 / (constants.R * constants.T_1); % kg/m3 section 1, page 6
constants.V_1 = constants.M_1 * sqrt(constants.gamma*constants.R*constants.T_1); % m/s, Section 1, page 10, or https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/mchkdrv.html

constants.m_dot = constants.rho_1 * constants.V_1 * constants.A_1; % kg / s

% calculate T_0 at position 1
constants.T_01 = (1 + ((constants.gamma-1)/2)*constants.M_1^2) * constants.T_1; % K, section 1, page 13
constants.T_02 = constants.T_01 + constants.Q_dot / (constants.m_dot * constants.c_p); % K, energy equation, section 2, page 5

% from problem statement, dT_0/dx is assumed constant from heating
constants.dT_0_dx = (constants.T_02 - constants.T_01) / constants.L; % K / m

dx = constants.dx;
L = constants.L;

xs = 0:dx:constants.L;
Ms_area     = zeros([1 length(xs)]);
Ms_friction = zeros([1 length(xs)]);
Ms_heat     = zeros([1 length(xs)]);
Ms_area(1)      = constants.M_1;
Ms_friction(1)  = constants.M_1;
Ms_heat(1)      = constants.M_1;

for i = 1:length(xs)-1
    x = xs(i);

    % calculate Area driven potential using runge kutta 4th order
    Ms_area(i+1) = myRungeKutta4(@dp_area_M,constants,Ms_area(i),x,dx);

    % calculate Friction driven potential using runge kutta 4th order
    Ms_friction(i+1) = myRungeKutta4(@dp_friction_M,constants,Ms_friction(i),x,dx);

    % calculate Heat driven potential using runge kutta 4th order
    Ms_heat(i+1) = myRungeKutta4(@dp_heat_M,constants,Ms_heat(i),x,dx);
end

figure
hold on 
plot(xs,Ms_area,"DisplayName",sprintf("Area Driven, M_L = %.5f",Ms_area(end)))
plot(xs,Ms_friction,"DisplayName",sprintf("Friction Driven, M_L = %.5f",Ms_friction(end)))
plot(xs,Ms_heat,"DisplayName",sprintf("Heat Driven, M_L = %.5f",Ms_heat(end)))
legend("Location","best")
xlabel("Duct Length [m]")
ylabel("Mach")
grid