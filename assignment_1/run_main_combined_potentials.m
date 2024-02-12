%% Notes
% CONVERGING DUCT
%
% Main run script for problem 1-3
%
% Note: change the dxs input for desired input for problem 2 or 3

%%

% cleanup
clear
clc
close all

% key inputs
constants = {};
constants.M_1 = 0.2;
%constants.dx = 0.25; % m

%dxs = [0.25,0.5,1.0]; % input for problem 2
dxs = [0.5]; % input for problem 3

line_styles = ["-","-","-"];

% initialize plots
f1 = figure(1);
ylabel('Mach')
xlabel('Duct Length [m]')
title(sprintf("Mach, M_1=%.1f",constants.M_1))
legend
grid

f2 = figure(2);
ylabel('Ratio')
xlabel('Duct Length [m]')
title(sprintf("T/T_1, M_1=%.1f",constants.M_1))
legend
grid

f3 = figure(3);
ylabel('Ratio')
xlabel('Duct Length [m]')
title(sprintf("P/P_1, M_1=%.1f",constants.M_1))
legend
grid

% outer loop to vary the dx step size
for j = 1:length(dxs)
    line_style = line_styles(j);
    constants.dx = dxs(j);
    disp(line_style)

    % constant inputs
    constants.gamma = 1.4;
    constants.alpha = deg2rad(1); % radians
    constants.R = 287; % J / kg * K
    constants.f_f = 0.025; % m
    constants.Q_dot = 500*1000; % J / s
    constants.L = 8.0; % m
    constants.c_p = (constants.gamma*constants.R) / (constants.gamma - 1); % section 1, page 8, cp constant for certain temperature range (ones we care about)
    
    % simulation flags
    constants.converging = true; % duct is converging
    constants.combined_potential = true; % multiple driving potentials
    
    % inlet conditions
    constants.r_1 = 0.5; % m
    constants.A_1 = pi*constants.r_1^2; % m^2
    constants.T_1 = 300; % K
    constants.P_1 = 100000; % Pa, kg/m*s2
    constants.rho_1 = constants.P_1 / (constants.R * constants.T_1); % kg/m3 section 1, page 6
    constants.V_1 = constants.M_1 * sqrt(constants.gamma*constants.R*constants.T_1); % m/s, Section 1, page 10, or https://www.grc.nasa.gov/www/k-12/VirtualAero/BottleRocket/airplane/mchkdrv.html
    constants.P_01 = constants.P_1 * ( 1+((constants.gamma-1)/2)*constants.M_1^2 )^(constants.gamma/(constants.gamma-1)); % section 2, page 7/13
    
    constants.m_dot1 = constants.rho_1 * constants.V_1 * constants.A_1; % kg / s
    
    % calculate T_0 at position 1
    constants.T_01 = (1 + ((constants.gamma-1)/2)*constants.M_1^2) * constants.T_1; % K, section 1, page 13
    constants.T_02 = constants.T_01 + constants.Q_dot / (constants.m_dot1 * constants.c_p); % K, energy equation, section 2, page 5
    
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
        
        M2_combined(i+1) = myRungeKutta4(@dp_combined_M,constants,M2_combined(i),x,dx);
        
    end
    
    M = sqrt(M2_combined);
    
    % calculate stagnation temperature across the duct using dT0_dx
    T_0x = constants.T_01 + constants.dT_0_dx * xs;

    % calculate temperature distribution
    num = 1 + ((constants.gamma-1)/2)*constants.M_1^2; % scalar
    den = 1 + ((constants.gamma-1)/2).*M2_combined; % array
    Tx_T1 = (T_0x/constants.T_01).*(num./den); % isentropic relation, section 2, page 32
    Tx = Tx_T1.*constants.T_1;
    
    % calculate stagnation temperature distribution
    T0x_T01 = T_0x / constants.T_01;
    
    % calculate velocity distribution
    Mx = sqrt(M2_combined);
    Vx_V1 = (Mx/constants.M_1).*sqrt(Tx_T1); % definition of Mach number, section 2, page 33
    V1_Vx = Vx_V1.^-1;
    Vx = Vx_V1.*constants.V_1;
    
    % calculate density distribution
    Ax = pi * (constants.r_1-xs*tan(constants.alpha)).^2; % for converging duct
    A1_Ax = constants.A_1 ./ Ax;
    rhox_rho1 = V1_Vx .* A1_Ax; % continuity equation, section 2, page 33
    
    % calculate pressure distribution
    Px_P1 = rhox_rho1 .* Tx_T1; % equation of state, section 2, page 33
    Px = Px_P1 * constants.P_1;
    
    % calculate stagnation pressure distribution
    % isentropic relation, section 2, page 33
    P_0x = Px .* ( 1+((constants.gamma-1)/2)*Mx.^2 ).^(constants.gamma/(constants.gamma-1));
    P0x_P01 = P_0x ./ constants.P_01;
    
    for i = 1:length(xs)
        fprintf('%8.2f%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n',xs(i),M(i),Vx_V1(i),Tx_T1(i),Px_P1(i),P0x_P01(i),T0x_T01(i))
    end
    
    figure(1)
    hold on 
    plot(xs,Mx,"DisplayName",sprintf("M_L = %.5f, dx=%.1f",Mx(end),constants.dx),'LineStyle',line_style)
    hold off
    
    figure(2)
    hold on 
    plot(xs,Tx_T1,"DisplayName",sprintf("T/T1 dx=%.1f",constants.dx),'LineStyle',line_style)
    hold off
    
    figure(3)
    hold on 
    plot(xs,Px_P1,"DisplayName",sprintf("P/P1 dx=%.1f",constants.dx),'LineStyle',line_style)
    hold off

    %figure(4)
    %hold on 
    %plot(xs,Vx_V1,"DisplayName",sprintf("V/V1"))
    %hold off

end