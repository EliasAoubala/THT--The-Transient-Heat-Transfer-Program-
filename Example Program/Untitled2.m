clear
clc
clf

% Script for Practice Utilisation of OneDimentionalHeatTransfer function
n_1=10; % Number of Points across the array 
n_2=20;

L_1=10*10^-3; % Thickness of Plate (mm)
L_2=20*10^-3; 

K_1=470; % Thermal Condusctivity (W/M*K) for mild steel.
K_2=19;

rho_1=2230; % Density of steel in kg/m^3
rho_2=8060;

Cp_1=717; % Specific Heat of steel (J/(kG K))
Cp_2=530;

Fo=1/4; % Fourier Number (tend to zero for higher accuracy)

t_end=60; % End Time (s)

T_int=293.15;   % Initial Temperature (K)
T_amb=T_int;    % Ambient Temperature (K)

%% Hot Gas Coefficient Calculations

% Bartz Hot Gas Calculation

%   - Constants Required for Bartz Coefficient

Pr=0.7536; % Prantl Number

Tc=1911.5; % Chamber Temperature (K)

gamma=1.3; % Gamma (N/D)

mue_g=0.000068759; % Dynamic Viscosity (Pa s)

Cp=1.2397*10^3; % Specific Heat (J/Kg*K)

Pc=10*10^5; % Chamber Pressure (Pa)

c_star=1100; % Throat Sonic Velocity (m/s)

At=6.63*10^-5; % Throat Area (m^2)

dt_throat=9.18*10^-3; % Throat Diameter (m)

%%   - Variable Parameters

M =0.1; % Mach Number

A_var = At*30; % Selected Area of investigation (m^2)

d_var= 2*(A_var/pi)^(1/2); % Selected Diameter of investigation (m)

hl=25; % 'Cold-Side' Heat Transfer Coefficient

CombP=[Pr, Tc, gamma, mue_g, Cp, Pc, c_star, At, dt_throat, M, A_var, d_var];

matArray=["Graphite", n_1, L_1, K_1, rho_1, Cp_1, 3800;
          "Stainless", n_2, L_2, K_2, rho_2, Cp_2, 1100];

BC=[Fo, t_end, hl, 10, 20, 30, 60];

OneDimentionalHeatTransfer(T_int, T_amb, CombP, matArray, BC)