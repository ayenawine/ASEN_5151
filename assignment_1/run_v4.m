%% Notes
% should have 6 total output plots
% CONVERGING DUCT

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
constants.converging = true;

% inlet conditions
constants.r_1 = 0.5; % m
constants.A_1 = pi*constants.r_1^2; % m^2
constants.M_1 = 0.2;
constants.T_1 = 300; % K
constants.P_1 = 100000; % Pa, kg/m*s2
constants.rho_1 = constants.P_1 / (constants.R * constants.T_1); % kg/m3 section 1, page 6
constants.V_1 = constants.M_1 * sqrt(constants.gamma*constants.R*constants.T_1); % m/s, Section 1, page 10, or https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/mchkdrv.html
constants.P_01 = constants.P_1 * ( 1+((constants.gamma-1)/2)*constants.M_1^2 )^(constants.gamma/(constants.gamma-1)); % section 2, page 7/13

constants.m_dot = constants.rho_1 * constants.V_1 * constants.A_1; % kg / s

% calculate T_0 at position 1
constants.T_01 = (1 + ((constants.gamma-1)/2)*constants.M_1^2) * constants.T_1; % K, section 1, page 13
constants.T_02 = constants.T_01 + constants.Q_dot / (constants.m_dot * constants.c_p); % K, energy equation, section 2, page 5

% from problem statement, dT_0/dx is assumed constant from heating
constants.dT_0_dx = (constants.T_02 - constants.T_01) / constants.L; % K / m

dx = constants.dx;
L = constants.L;

xs = 0:dx:constants.L;

M2_area = zeros([1 length(xs)]);
M2_fanno = zeros([1 length(xs)]);
M2_rayleigh = zeros([1 length(xs)]);
M2_combined = zeros([1 length(xs)]);

M2_area(1) = constants.M_1^2;
M2_fanno(1) = constants.M_1^2;
M2_rayleigh(1) = constants.M_1^2;
M2_combined(1) = constants.M_1^2;

for i = 1:length(xs)-1
    x = xs(i);

    % NOTE TO SELF:
    % ODE is M2 with respect to dx, so M2 has to be input into the runge
    % kutta, not M

    % calculate Area driven potential using runge kutta 4th order
    %M2_area(i+1) = myRungeKutta4(@dp_area_M,constants,M2_area(i),x,dx);

    % calculate Fanno driven potential using runge kutta 4th order
    %M2_fanno(i+1) = myRungeKutta4(@dp_friction_M,constants,M2_fanno(i),x,dx);

    % calculate Rayleigh driven potential using runge kutta 4th order
    %M2_rayleigh(i+1) = myRungeKutta4(@dp_heat_M,constants,M2_rayleigh(i),x,dx);

    % Calculate potential from all 3 flows
    %inputs = {constants,M2_rayleigh(i),x,dx};
    %dp_combined_M = @(inputs) dp_area_M(inputs)+dp_friction_M(inputs)+dp_heat_M(inputs);
    
    M2_combined(i+1) = myRungeKutta4(@dp_combined_M,constants,M2_combined(i),x,dx);
    
end

% calculate temperature distribution (section 2, page 32)
T_0x = constants.T_01 + constants.dT_0_dx * xs;
num = 1 + ((constants.gamma-1)/2)*constants.M_1^2; % scalar
den = 1 + ((constants.gamma-1)/2).*M2_combined; % array
Tx_T1 = (T_0x/constants.T_01).*(num./den);
Tx = Tx_T1.*constants.T_1;

% calculate velocity distribution
Mx = sqrt(M2_combined);
Vx_V1 = (Mx./constants.M_1).*sqrt(Tx/constants.T_1);
Vx = Vx_V1.*constants.V_1;

% calculate density distribution
Ax = pi.*(constants.r_1-xs*tan(constants.alpha)).^2; % for converging duct
A1_Ax = constants.A_1 ./ Ax;
rhox_rho1 = A1_Ax ./ Vx_V1;

% calculate pressure distribution
Px_P1 = rhox_rho1 ./ Tx_T1;

% calculate pressure distribution (section 2, page 33)
%P_0x = constants.P_1 * ( 1+((constants.gamma-1)/2)*constants.M_1^2 )^(constants.gamma/(constants.gamma-1));

figure
hold on 
plot(xs,sqrt(M2_combined),"DisplayName",sprintf("Combined Potential M_L = %.5f",sqrt(M2_combined(end))))
legend
grid

figure
hold on 
plot(xs,Tx_T1,"DisplayName",sprintf("T/T1"))
legend
grid

figure
hold on 
plot(xs,Vx_V1,"DisplayName",sprintf("V/V1"))
legend
grid