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

%% Aircraft Data (Learjet 24 Aircraft)
% Geometric data
S = .55;                            %Wing surface area [m^2]
c_bar = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]

% Flight Conditions Data
g = 9.8;
altt = 100;                       %Altitude [m]
Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]
%alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

%% Trim

V_a0 = 30;
alpha0 = 0;
delta_e0 = deg2rad(4);
delta_t0 = 0;

w0 = 0;
q0 = 0;


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
Gamma_6 = J_xz/Jy;
Gamma_7 = (J_x*(J_x-J_y) + J_xz^2)/Gamma;
Gamma_8 = J_x/Gamma;

% Steady level trim
u0 = V;                             %Initial velocity [ft/sec]
theta0 = deg2rad(0);                %Initial pitch angle [rad]

%% Prop properties
S_prop = 0.2027;
k_motor = 80;
C_prop = 1;

%% Longitudinal Aerodynamic Coefficients
% Steady State
C_L1 = 0.41;
C_D1 = 0.0335;
C_m1 = 0;
C_Tx1 = 0.0335;
C_mT1 = 0;

% Stability Derivatives
C_D0 = 0.03;
C_D_alpha = 0.3;
C_D_q = 0;
C_L0 = 0.28;
C_L_alpha = 3.45;
C_L_q = 0;
C_m0 = -0.02338;
C_mu = 0;%0.05;?
C_m_alpha = -0.38;

C_m_q = -3.6;



% Control Derivatives
C_D_delta_e = 0;
C_Dih = 0;%?
C_L_delta_e = -0.36;
C_Lih = 0;%0.94;%?
C_m_delta_e = -0.5;
C_mih = 0;% -2.5;?

%% Lateral Directional Aerodynamic Coefficients
% Stability Derivatives
C_lbeta = -0.12;
C_lp = -0.26;
C_lr = 0.14;
C_Ybeta = -0.98;
C_Yp = 0;
C_Yr = 0;
C_nbeta = 0.25;
C_nTbeta = 0;%?
C_np = 0.022;
C_nr = -0.35;

% Control Derivatives
C_ldeltaa = 0.08;
C_ldeltar = 0.105;
C_Ydeltaa = 0;
C_Ydeltar = -0.17;
C_ndeltaa = 0.06;
C_ndeltar = 0;%-0.074;?
%% Standard Atmospheric Calculations

[T, a, P, rho] = atmosisa(altt);
    
%% Nondimensional Stability Derivatives
%From "translation" key in HW3
%C_xu = -(C_Du+2*C_D1);
C_X_0 = -C_D0*cos(alpha0) + C_L0*sin(alpha0);
C_X_alpha = -C_D_alpha*cos(alpha0) + C_L_alpha*sin(alpha0);
%C_xalpha_dot = 0;
C_X_q = -C_D_q*cos(alpha0) + C_L_q*sin(alpha0);
C_X_delta_e = -C_D_delta_e*cos(alpha0) + C_L_delta_e*sin(alpha0);
%C_X_deltat = C_Txu+2*C_Tx1;

% C_Z_u = -(C_Lu+2*C_L1);
% C_Z_alpha = -(C_Lalpha+C_D1);
% C_zalpha_dot = -C_Lalpha_dot;
% C_Z_q = -C_Lq;
% C_Z_delta_e = -C_Ldeltae;
% C_zdeltap = 0;

C_Z_0 = -C_D0*sin(alpha0) + C_L0*cos(alpha0);
C_Z_alpha = -C_D_alpha*sin(alpha0) + C_L_alpha*cos(alpha0);
C_Z_q = -C_D_q*sin(alpha0) + C_L_q*cos(alpha0);
C_Z_delta_e = -C_D_delta_e*sin(alpha0) + C_L_delta_e*cos(alpha0);

%C_mu = C_mu+2*C_m1;
%C_m_alpha = C_malpha;
%C_malpha_dot = C_malpha_dot;
%C_m_q = C_m_q;
%C_m_delta_e = C_m_delta_e;
%C_mdeltap = C_mTu+2*C_mT1;

% Nondimensional
%C_w0 = w/(1/2*rho*u0^2*S);          %Non-dimensional weight [unitless]

% Longitudinal Dimensional Derivatives
Xu = (u0*rho*S)/m * (C_X_0 + C_X_alpha*alpha0 + C_X_delta_e*delta_e0) - (rho*S*w0*C_X_alpha)/2*m + (rho*S*c_bar*C_X_q*u0*q0)/4*m*V_a0 - (rho*S_prop*C_prop*u0)/m;
Xw = -q0 + (w0*rho*S)/m * (C_X_0 + C_X_alpha*alpha0 + C_X_delta_e*delta_e0) + (rho*S*c_bar*C_X_q*w0*q0)/4*m*V_a0 + (rho*S*u0*C_X_alpha)/2*m - (rho*S_prop*C_prop*w0)/m;
Xq = -w0 + (rho*V_a0*c_bar*S*C_X_q)/4*m;
%Xw_dot = 1/4*rho*c_bar*S*C_xalpha_dot;

Zu = q0 + (u0*rho*S)/m * (C_Z_0 + C_Z_alpha*alpha0 + C_Z_delta_e*delta_e0) - (rho*S*w0*C_Z_alpha)/2*m + (rho*S*c_bar*C_Z_q*u0*q0)/4*m*V_a0;
Zw = (w0*rho*S)/m * (C_Z_0 + C_Z_alpha*alpha0 + C_Z_delta_e*delta_e0) - (rho*S*u0*C_Z_alpha)/2*m + (rho*S*c_bar*C_Z_q*w0*q0)/4*m*V_a0;
Zq = u0 + (rho*V_a0*S*C_Z_q*q0)/4*m;
%Zw_dot = 1/4*rho*c_bar*S*C_zalpha_dot;

Mu = ((u0*rho*S*c_bar)/J_y)*(C_m0 + C_m_alpha*alpha0 + C_m_delta_e*delta_e0) - (rho*S*c_bar*C_m_alpha*u0)/2*J_y + (rho*S*c_bar^2 *q0*u0)/(4*J_y*V_a0);
Mw = ((w0*rho*S*c_bar)/J_y)*(C_m0 + C_m_alpha*alpha0 + C_m_delta_e*delta_e0) + (rho*S*c_bar*C_m_alpha*u0)/2*J_y + (rho*S*c_bar^2 *q0*u0)/(4*J_y*V_a0);
Mq = (1/4*rho*V_a0*c_bar^2*S*C_m_q)/J_y;
%Mw_dot = 1/4*rho*c_bar^2*S*C_malpha_dot;

%Lateral Dimensional Derivatives
Yv = 1/2*rho*u0*S*C_Ybeta;
Yp = 1/4*rho*u0*b*S*C_Yp;
Yr = 1/4*rho*u0*b*S*C_Yr;

Lv = 1/2*rho*u0*b*S*C_lbeta;
Lp = 1/4*rho*u0*b^2*S*C_lp;
Lr = 1/4*rho*u0*b^2*S*C_lr;

Nv = 1/2*rho*u0*b*S*C_nbeta;
Np = 1/4*rho*u0*b^2*S*C_np;
Nr = 1/4*rho*u0*b^2*S*C_nr;

%% Dimensional Control Derivatives
%Longitudinal Dimensional Control Derivatives
X_deltae = (C_X_delta_e*1/2*rho*V_a^2*S)/m;
X_deltat = (rho*S_prop*C_prop*k_motor^2*delta_t0)/m;

Z_deltae = (C_Z_delta_e*1/2*rho*V_a0^2*S)/m;
%Z_deltap = C_zdeltap*1/2*rho*u0^2*S;

M_deltae = (C_m_delta_e*1/2*rho*V_a0^2*S*c_bar)/J_y;
%M_deltap = C_mdeltap*1/2*rho*u0^2*S*c_bar;

%Lateral Dimensional Control Derivatives
Y_deltaa = C_Ydeltaa*1/2*rho*u0^2*S;
Y_deltar = C_Ydeltar*1/2*rho*u0^2*S;

L_deltaa = C_ldeltaa*1/2*rho*u0^2*S*b;
L_deltar = C_ldeltar*1/2*rho*u0^2*S*b;

N_deltaa = C_ndeltaa*1/2*rho*u0^2*S*b;
N_deltar = C_ndeltar*1/2*rho*u0^2*S*b;

%% Longitudinal Linear Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [u,w,q,theta,h]'
Long1A = Xu;                                            %Row 1, Column 1
Long2A = Xw;                                            %Row 1, Column 2
Long3A = Xq;                                               %Row 1, Column 3
Long4A = -g*cos(theta0);                                  %Row 1, Column 4
Long1_5A = 0;
Long5A = Zu;                                   %Row 2, Column 1
Long6A = Zw;                                   %Row 2, Column 2
Long7A = Zq;                          %Row 2, Column 3
Long8A = -g*sin(theta0);                   %Row 2, Column 4
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
Long5_1A = sin(theta0);
Long5_2A = -cos(theta0);
Long5_3A = 0;
Long5_4A = u0*cos(theta0)+w0*sin(theta0);
Long5_5A = 0;


A_long = [Long1A, Long2A, Long3A, Long4A Long1_5A;
          Long5A, Long6A, Long7A, Long8A Long2_5A;
          Long9A, Long10A, Long11A, Long12A Long3_5A;
          Long13A, Long14A, Long15A, Long16A Long4_5A;
          Long5_1A Long5_2A Long5_3A Long5_4A Long5_5A]
      
B_long = [X_deltae X_deltat;
          Z_deltae 0;
          M_deltae 0;
          0 0;
          0 0;]
      
%% Lateral Linear Model
Lat1A = Y_v;
Lat2A = Y_p;
Lat3A = Y_r;
Lat4A = g*cos(theta0)*cos(phi0);
Lat1_5A = 0;
Lat5A = L_v;
Lat6A = L_p;
Lat7A = L_r;
Lat8A = 0;
Lat2_5A = 0;
Lat9A = N_v;
Lat10A = N_p;
Lat11A = N_r;
Lat12A = 0;
Lat3_5A = 0;
Lat13A = 0;
Lat14A = 1;
Lat15A = cos(phi0)*tan(theta0);
Lat16A = q0*cos(phi0)*tan(theta0) - r0*sin(phi0)*tan(theta0);
Lat4_5A = 0;
Lat5_1A = 0;
Lat5_2A = 0;
Lat5_3A = cos(phi0)*1/cos(theta0);
Lat5_4A = p0*cos(phi0)*1/cos(theta0) - r0*sin(phi0)*1/cos(theta0);
Lat5_5A = 0;

A_lat = [Lat1A, Lat2A, Lat3A, Lat4A Lat1_5A;
          Lat5A, Lat6A, Lat7A, Lat8A Lat2_5A;
          Lat9A, Lat10A, Lat11A, Lat12A Lat3_5A;
          Lat13A, Lat14A, Lat15A, Lat16A Lat4_5A;
          Lat5_1A Lat5_2A Lat5_3A Lat5_4A Lat5_5A]
      
B_lat = [Y_delta_a                              Y_delta_r;
         L_delta_a  L_delta_r;
         N_delta_a  (N_deltar/Iz_prime)+Izx_prime*L_deltar;
         0                                       0]

     
%f_p_x = .5*rho*S_prop*C_prop*((k_motor*del_t)^2 - V_a^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
% del_a
% del_p
% del_r_r
% del_r_l

% [el_rud] = [1 1; -1 1]*[del_r_r; del_r_l];
% 
% del_e = el_rud(1);
% del_r = el_rud(2);
% 
% X_long = [p_x;
%           p_z;
%           u_t;
%           w_t;
%           theta;
%           q];
% 
%       
% X_lat = [p_y;
%          v_t;
%          psi;
%          phi;
%          p;
%          r];
%      
% U_long = [del_e;
%           del_p];
% 
% U_lat = [del_a;
%          del_r];
%      
% X_dot_long = A_long*X_long + B_long*U_long
% 
% X_dot_lat = A_lat*X_lat + B_lat*U_lat
% 
% 
% 
% f = [-m*g*sin(theta);
%      m*g*cos(theta)*sin(phi);
%      m*g*cos(theta)*cos(phi);]
%      
% 
