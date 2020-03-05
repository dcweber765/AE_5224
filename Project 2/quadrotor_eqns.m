clc; clear all;
%% Hummingbird Parameters
I_xx = 2.32e-3; %kg/m^2
I_yy = 2.32e-3; %kg/m^2
% I_xx = I_yy -> Rotor craft is symmetric
I_zz = 4.00e-3; %kg/m^2
m = 0.5;        %kg
L = 0.175;      %m 
k_M = 1.5e-9;   %Nm/rpm^2
k_F = 6.11e-8;  %N/rpm^2
k_m = 20;       %s^-1
gamma=k_M/k_F;
g = 9.81;       %m/s^2
W = m*g;        %N
%% Flight conditions
altt = 100;
[T, a, P, rho] = atmosisa(altt);
%% Rotor Forces
F1 = k_F*omega_1^2;
F2 = k_F*omega_2^2;
F3 = k_F*omega_3^2;
F4 = k_F*omega_4^2;
%% Inertial Matrix
I = [I_xx 0 0;0 I_yy 0;0 0 I_zz];
%% Position rates
R_v_b = [cos(theta)*cos(psi),sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi),cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
           cos(theta)*sin(psi),sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi),cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
           sin(theta)         ,-sin(phi)*cos(theta)                          ,-cos(phi)*cos(theta)                         ];
V_b = [u v w]';
p_dot = R_v_b_T*V_b;
p_x_dot = p_dot(1);
p_y_dot = p_dot(2);
p_z_dot = p_dot(3);
%% Wind rates
u1 = F1+F2+F3+F4;
% SCD0 =  w*0.95; %Piazza -> S*CD0 should be less than 10% of the weight
% Drag = (1/2)*rho*V^2*SCD0;
V_dot = [0 0 -g]'+ (1/m)*inv(R_v_b)*[0 0 u1];
u_dot = V_dot(1);
v_dot = V_dot(2);
w_dot = V_dot(3);
%% Euler rates
omega_b_b = [p q r]';
m1 = [cos(theta) 0 -cos(phi)*sin(theta);0 1 sin(phi);sin(theta) 0 cos(phi)*cos(theta)];
euler_dot = inv(m1)*omega_b_b;
phi_dot = euler_dot(1);
theta_dot = euler_dot(2);
psi_dot = euler_dot(3);
%% Angular rates
u2 = [0 L 0 -L; -L 0 L 0; gamma -gamma gamma -gamma]*[F1 F2 F3 F4]';
ang_dot = inv(I)* u2-cross(omega_b_b,I*omega_b_b);
p_dot = ang_dot(1);
q_dot = ang_dot(2);
r_dot = ang_dot(3);