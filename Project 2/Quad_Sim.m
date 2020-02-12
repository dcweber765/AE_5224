clear all; close all;

%% Quad Rotor Sim

%% Flight Conditions Data

g = 9.8;
altt = 100;                       %Altitude [m]
V = 30;                            %True airspeed [m/sec]
[T, SoS, P, rho] = atmosisa(altt); %Standard Atmospheric Calculations

%% Aircraft Data

J_x = 2.32E-3; %kg*m^-2
J_y = 2.32E-3; 
J_z = 4E-3;

m = .5;%kg
L = .175;%m

psi_0 = 0;

a = [g*cos(psi_0) g*sin(psi_0) 0;...
    g*sin(psi_0) -g*cos(psi_0) 0;...
    0 0 0];

A = [zeros(3,3), eye(3),zeros(3,6);...
    zeros(3,6), a, zeros(3,3);...
    zeros(3,3), zeros(3,6), eye(3);...
    zeros(3,3) zeros(3,6), zeros(3,3)];


B = [zeros(5,4);...
    1/m 0 0 0;...
    zeros(3,4);...
    0 1/J_x 0 0;...
    0 0 1/J_y 0;...
    0 0 0 1/J_z];
