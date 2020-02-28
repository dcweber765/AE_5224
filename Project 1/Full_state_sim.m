%% FROM ADC - B18 - HW3

clear all; clear all;

%% STATE
      p_x = 0;
      p_y = 0;
      p_z = 0;
      u_t = 0;
      v_t = 0;
      w_t = 0;
      u_b = 0;
      v_b = 0;
      w_b = 0;
      V_a = 0;
      alpha = 0;
      beta = 0;
      psi = 0;
      theta = 0;
      phi = 0;
      p = 0;
      q = 0;
      r = 0;
      
p_t = [p_x;p_y;p_z];
V_t = [u_t;v_t;w_t];
V_b = [u_b;v_b;w_b];
euler_ang = [psi;theta;phi];
omega_tb_b = [p;q;r];

X = [p_t;
    V_t;
    V_b;
    V_a;
    alpha;
    beta;
    euler_ang;
    omega_tb_b;]

%% Aircraft Data (UAV)
% Geometric data
S = .55;                            %Wing surface area [m^2]
c_bar = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]
%h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]


%% Flight Conditions Data

g = 9.8;
altt = 100;                       %Altitude [m]
Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
[T, a, P, rho] = atmosisa(altt); %Standard Atmospheric Calculations
%alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

%% Trim

V_a_star = 30;
alpha_star = 0.0382;
beta_star = 0;
delta_e_star = -.0758;
delta_t_star = .3955;
delta_a_star = 0;
delta_r_star = 0;

u_star = 29.9781;                             %Initial velocity [ft/sec]
v_star = 0;
w_star = 1.1451;

p_star = .0382;
q_star = 0;
r_star = 0;

psi_star = 0;
theta_star = deg2rad(0);                %Initial pitch angle [rad]
phi_star = 0;

% Mass and Inertial Data

m = 13.5;                            %Mass of the aircraft [kg]
w = m*g;                          %Weight of the aircraft [N]
J_x = 0.8244;                         %Moment of inertia x-axis [kg*m^2]
J_y = 1.135;                         %Moment of inertia y-axis [kg*m^2]
J_z = 1.759;                         %Moment of inertia z-axis [kg*m^2]
J_zx = 0.1204;                         %Moment of inertia zx-axis [kg*m^2]
J_xz = J_zx;                          %Moment of inertia xz-axis [kg*m^2]

Gamma = J_x*J_z - J_xz^2;

Gamma_1 = J_xz*(J_x-J_y+J_z)/Gamma;
Gamma_2 = (J_z*(J_z-J_y) + J_xz^2)/Gamma;
Gamma_3 = J_z/Gamma;
Gamma_4 = J_xz/Gamma;
Gamma_5 = (J_z - J_x)/J_y;
Gamma_6 = J_xz/J_y;
Gamma_7 = (J_x*(J_x-J_y) + J_xz^2)/Gamma;
Gamma_8 = J_x/Gamma;

% Steady level trim

%% Prop properties
S_prop = 0.2027;
k_motor = 80;
C_prop = 1;

%% Longitudinal Aerodynamic Coefficients
% Steady State
% C_L1 = 0.41;
% C_D1 = 0.0335;
% C_m1 = 0;
% C_Tx1 = 0.0335;
% C_mT1 = 0;

% Stability Derivatives
C_D_0 = 0.03;
C_D_alpha = 0.3;
C_D_q = 0;
C_L0 = 0.28;
C_L_alpha = 3.45;
C_L_q = 0;
C_m_0 = -0.02338;
C_m_alpha = -0.38;
C_m_q = -3.6;

C_D_p = 0.0437;

% Control Derivatives
C_D_delta_e = 0;
C_L_delta_e = -0.36;
C_m_delta_e = -0.5;

%% Lateral Directional Aerodynamic Coefficients
% Stability Derivatives
C_l_0 = 0;
C_l_beta = -0.12;
C_l_p = -0.26;
C_l_r = 0.14;
C_Y_0 = 0;
C_Y_beta = -0.98;
C_Y_p = 0;
C_Y_r = 0;
C_n_0 = 0;
C_n_beta = 0.25;
C_n_p = 0.022;
C_n_r = -0.35;

% Control Derivatives
C_l_delta_a = 0.08;
C_l_delta_r = 0.105;
C_Y_delta_a = 0;
C_Y_delta_r = -0.17;
C_n_delta_a = 0.06;
C_n_delta_r = -0.032;

    
%% Nondimensional Stability Derivatives

C_X_0 = -C_D_0*cos(alpha_star) + C_L0*sin(alpha_star);
C_X_alpha = -C_D_alpha*cos(alpha_star) + C_L_alpha*sin(alpha_star);
C_X_q = -C_D_q*cos(alpha_star) + C_L_q*sin(alpha_star);
C_X_delta_e = -C_D_delta_e*cos(alpha_star) + C_L_delta_e*sin(alpha_star);

C_Z_0 = -C_D_0*sin(alpha_star) + C_L0*cos(alpha_star);
C_Z_alpha = -C_D_alpha*sin(alpha_star) + C_L_alpha*cos(alpha_star);
C_Z_q = -C_D_q*sin(alpha_star) + C_L_q*cos(alpha_star);
C_Z_delta_e = -C_D_delta_e*sin(alpha_star) + C_L_delta_e*cos(alpha_star);

C_p_0 =       Gamma_3*C_l_0 + Gamma_4*C_n_0;
C_p_beta =    Gamma_3*C_l_beta + Gamma_4*C_n_beta;
C_p_p =       Gamma_3*C_l_p + Gamma_4*C_n_p;
C_p_r =       Gamma_3*C_l_r + Gamma_4*C_n_r;
C_p_delta_a = Gamma_3*C_l_delta_a + Gamma_4*C_n_delta_a;
C_p_delta_r = Gamma_3*C_l_delta_r + Gamma_4*C_n_delta_r;

C_r_0 =       Gamma_4*C_l_0 + Gamma_8*C_n_0;
C_r_beta =    Gamma_4*C_l_beta + Gamma_8*C_n_beta;
C_r_p =       Gamma_4*C_l_p + Gamma_8*C_n_p;
C_r_r =       Gamma_4*C_l_r + Gamma_8*C_n_r;
C_r_delta_a = Gamma_4*C_l_delta_a + Gamma_8*C_n_delta_a;
C_r_delta_r = Gamma_4*C_l_delta_r + Gamma_8*C_n_delta_r;

% Longitudinal Dimensional Derivatives
Xu = (u_star*rho*S)/m * (C_X_0 + C_X_alpha*alpha_star + C_X_delta_e*delta_e_star) - (rho*S*w_star*C_X_alpha)/2*m + (rho*S*c_bar*C_X_q*u_star*q_star)/4*m*V_a_star - (rho*S_prop*C_prop*u_star)/m;
Xw = -q_star + (w_star*rho*S)/m * (C_X_0 + C_X_alpha*alpha_star + C_X_delta_e*delta_e_star) + (rho*S*c_bar*C_X_q*w_star*q_star)/4*m*V_a_star + (rho*S*u_star*C_X_alpha)/2*m - (rho*S_prop*C_prop*w_star)/m;
Xq = -w_star + (rho*V_a_star*c_bar*S*C_X_q)/4*m;

Zu = q_star + (u_star*rho*S)/m * (C_Z_0 + C_Z_alpha*alpha_star + C_Z_delta_e*delta_e_star) - (rho*S*w_star*C_Z_alpha)/2*m + (rho*S*c_bar*C_Z_q*u_star*q_star)/4*m*V_a_star;
Zw = (w_star*rho*S)/m * (C_Z_0 + C_Z_alpha*alpha_star + C_Z_delta_e*delta_e_star) - (rho*S*u_star*C_Z_alpha)/2*m + (rho*S*c_bar*C_Z_q*w_star*q_star)/4*m*V_a_star;
Zq = u_star + (rho*V_a_star*S*C_Z_q*q_star)/4*m;

Mu = ((u_star*rho*S*c_bar)/J_y)*(C_m_0 + C_m_alpha*alpha_star + C_m_delta_e*delta_e_star) - (rho*S*c_bar*C_m_alpha*u_star)/2*J_y + (rho*S*c_bar^2 *q_star*u_star)/(4*J_y*V_a_star);
Mw = ((w_star*rho*S*c_bar)/J_y)*(C_m_0 + C_m_alpha*alpha_star + C_m_delta_e*delta_e_star) + (rho*S*c_bar*C_m_alpha*u_star)/2*J_y + (rho*S*c_bar^2 *q_star*u_star)/(4*J_y*V_a_star);
Mq = (1/4*rho*V_a_star*c_bar^2*S*C_m_q)/J_y;

%Lateral Dimensional Derivatives
Yv = ((rho*S*b*v_star)/(4*m*V_a_star))*(C_Y_p*p_star + C_Y_r*r_star) + ((rho*S*v_star)/m)*(C_Y_0 + C_Y_beta*beta_star + C_Y_delta_a*delta_a_star + C_Y_delta_r*delta_r_star) + ((rho*S*C_Y_beta)/(2*m))*sqrt(u_star^2+w_star^2);
Yp = w_star + (rho*V_a_star*S*b)/(4*m)*C_Y_p;
Yr = u_star + (rho*V_a_star*S*b)/(4*m)*C_Y_r;

Lv = ((rho*S*b*v_star)/(4*V_a_star))*(C_p_p*p_star + C_p_r*r_star) + (rho*S*v_star)*(C_p_0 + C_p_beta*beta_star + C_p_delta_a*delta_a_star + C_p_delta_r*delta_r_star) + ((rho*S*C_p_beta)/(2))*sqrt(u_star^2+w_star^2);
Lp = Gamma_1*q_star+ (rho*V_a_star*S*b^2)/4 * C_p_p;
Lr = -Gamma_2*q_star+ (rho*V_a_star*S*b^2)/4 * C_p_r;

Nv = ((rho*S*b*v_star)/(4*V_a_star))*(C_r_p*p_star + C_r_r*r_star) + (rho*S*v_star)*(C_r_0 + C_r_beta*beta_star + C_r_delta_a*delta_a_star + C_r_delta_r*delta_r_star) + ((rho*S*C_r_beta)/2)*sqrt(u_star^2+w_star^2);
Np = Gamma_7*q_star+ (rho*V_a_star*S*b^2)/4 * C_r_p;
Nr = -Gamma_1*q_star+ (rho*V_a_star*S*b^2)/4 * C_r_r;

%% Dimensional Control Derivatives
%Longitudinal Dimensional Control Derivatives
X_delta_e = (C_X_delta_e*1/2*rho*V_a^2*S)/m;
X_delta_t = (rho*S_prop*C_prop*k_motor^2*delta_t_star)/m;

Z_delta_e = (C_Z_delta_e*1/2*rho*V_a_star^2*S)/m;

M_delta_e = (C_m_delta_e*1/2*rho*V_a_star^2*S*c_bar)/J_y;

%Lateral Dimensional Control Derivatives
Y_delta_a = (C_Y_delta_a*1/2*rho*V_a_star^2*S)/m;
Y_delta_r = (C_Y_delta_r*1/2*rho*V_a_star^2*S)/m;

L_delta_a = C_p_delta_a*1/2*rho*V_a_star^2*S*b;
L_delta_r = C_p_delta_r*1/2*rho*V_a_star^2*S*b;

N_delta_a = C_r_delta_a*1/2*rho*V_a_star^2*S*b;
N_delta_r = C_r_delta_r*1/2*rho*V_a_star^2*S*b;

%% Longitudinal Linear Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u,w,q,theta,h]'
Long1A = Xu;                                            %Row 1, Column 1
Long2A = Xw;                                            %Row 1, Column 2
Long3A = Xq;                                               %Row 1, Column 3
Long4A = -g*cos(theta_star);                                  %Row 1, Column 4
Long1_5A = 0;
Long5A = Zu;                                   %Row 2, Column 1
Long6A = Zw;                                   %Row 2, Column 2
Long7A = Zq;                          %Row 2, Column 3
Long8A = -g*sin(theta_star);                   %Row 2, Column 4
Long2_5A = 0;
Long9A = Mu;            %Row 3, Column 1
Long10A = Mw;           %Row 3, Column 2
Long11A = Mq;    %Row 3, Column 3
Long12A = 0;      %Row 3, Column 4
Long3_5A = 0;
Long13A = 0;                                              %Row 4, Column 1
Long14A = 0;                                              %Row 4, Column 2
Long15A = 1;                                              %Row 4, Column 3
Long16A = 0;                                              %Row 4, Column 4
Long4_5A = 0;
Long5_1A = sin(theta_star);
Long5_2A = -cos(theta_star);
Long5_3A = 0;
Long5_4A = u_star*cos(theta_star)+w_star*sin(theta_star);
Long5_5A = 0;


A_long = [Long1A, Long2A, Long3A, Long4A Long1_5A;
          Long5A, Long6A, Long7A, Long8A Long2_5A;
          Long9A, Long10A, Long11A, Long12A Long3_5A;
          Long13A, Long14A, Long15A, Long16A Long4_5A;
          Long5_1A Long5_2A Long5_3A Long5_4A Long5_5A]
      
B_long = [X_delta_e X_delta_t;
          Z_delta_e 0;
          M_delta_e 0;
          0 0;
          0 0;]
      
%% Lateral Linear Model
Lat1A = Yv;
Lat2A = Yp;
Lat3A = Yr;
Lat4A = g*cos(theta_star)*cos(phi_star);
Lat1_5A = 0;
Lat5A = Lv;
Lat6A = Lp;
Lat7A = Lr;
Lat8A = 0;
Lat2_5A = 0;
Lat9A = Nv;
Lat10A = Np;
Lat11A = Nr;
Lat12A = 0;
Lat3_5A = 0;
Lat13A = 0;
Lat14A = 1;
Lat15A = cos(phi_star)*tan(theta_star);
Lat16A = q_star*cos(phi_star)*tan(theta_star) - r_star*sin(phi_star)*tan(theta_star);
Lat4_5A = 0;
Lat5_1A = 0;
Lat5_2A = 0;
Lat5_3A = cos(phi_star)*1/cos(theta_star);
Lat5_4A = p_star*cos(phi_star)*1/cos(theta_star) - r_star*sin(phi_star)*1/cos(theta_star);
Lat5_5A = 0;

A_lat = [Lat1A, Lat2A, Lat3A, Lat4A Lat1_5A;
          Lat5A, Lat6A, Lat7A, Lat8A Lat2_5A;
          Lat9A, Lat10A, Lat11A, Lat12A Lat3_5A;
          Lat13A, Lat14A, Lat15A, Lat16A Lat4_5A;
          Lat5_1A Lat5_2A Lat5_3A Lat5_4A Lat5_5A];
      
B_lat = [Y_delta_a Y_delta_r;
         L_delta_a  L_delta_r;
         N_delta_a  N_delta_r;
         0 0;
         0 0];

     
%f_p_x = .5*rho*S_prop*C_prop*((k_motor*del_t)^2 - V_a^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
% Delta_del_a
% Delta_del_p
% Delta_del_r_r
% Delta_del_r_l

[el_rud] = [1 1; -1 1]*[del_r_r; del_r_l];

Delta_del_e = el_rud(1);
Delta_del_r = el_rud(2);

Delta_X_long = [Delta_u_t; Delta_w_t; Delta_q; Delta_theta; Delta_p_z];

Delta_X_lat = [Delta_v_t; Delta_p; Delta_r; Delta_phi; Delta_psi];
     
Delta_U_long = [Delta_delta_e; Delta_delta_p];

Delta_U_lat = [Delta_delta_a; Delta_delta_r];
     
Delta_X_dot_long = A_long*Delta_X_long + B_long*Delta_U_long

Delta_X_dot_lat = A_lat*Delta_X_lat + B_lat*Delta_U_lat

Delta_X_long_plus = Delta_X_dot_long*dT
Delta_X_lat_plus = Delta_X_dot_lat*dT

Delta_X = [Delta_X_long_plus(1)*dT;...
           Delta_X_lat_plus(1)*dT;...
           Delta_X_long_plus(5);...
           Delta_X_long_plus(1);...
           Delta_X_lat_plus(1);...
           Delta_X_long_plus(2);...];
           
           
Long_sys = ss(A_long,B_long,eye(5),0);
Lat_sys = ss(A_lat,B_lat,zeros(5,5),0);
x0 = [0;0;0;0;0]
t = linspace(0,1,100);
u  = [delta_e_star*ones(1,100); delta_t_star*ones(1,100)];
lsim(Long_sys, u, t, x0)

% 
% f = [-m*g*sin(theta);
%      m*g*cos(theta)*sin(phi);
%      m*g*cos(theta)*cos(phi);]
%      
% 
