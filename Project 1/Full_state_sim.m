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
    V_a;X
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
%Mach = 0.7;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]
%alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

% Mass and Inertial Data

m = 13.5;                            %Mass of the aircraft [kg]
w = m*g;                          %Weight of the aircraft [N]
Jx = 0.8244;                         %Moment of inertia x-axis [kg*m^2]
Jy = 1.135;                         %Moment of inertia y-axis [kg*m^2]
Jz = 1.759;                         %Moment of inertia z-axis [kg*m^2]
Jzx = 0.1204;                         %Moment of inertia zx-axis [kg*m^2]
Jxz = Jzx;                          %Moment of inertia xz-axis [kg*m^2]

Ix_prime = (Jx*Jz-Jzx^2)/Jz;
Iz_prime = (Jx*Jz-Jzx^2)/Jx;
Izx_prime = Jzx/(Jx*Jz-Jxz^2);

% Steady level trim
u0 = V;                             %Initial velocity [ft/sec]
theta0 = deg2rad(0);                %Initial pitch angle [rad]
%% Longitudinal Aerodynamic Coefficients
% Steady State
C_L1 = 0.41;
C_D1 = 0.0335;
C_m1 = 0;
C_Tx1 = 0.0335;
C_mT1 = 0;

% Stability Derivatives
C_D0 = 0.03;
C_Du = 0;%0.104;?
C_Dalpha = 0.3;
C_Txu = 0;%-0.07;?
C_L0 = 0.28;
C_Lu = 0;%0.4;?
C_Lalpha = 3.45;
C_Lalpha_dot = 0;%2.2;?
C_Lq = 0;
C_m0 = -0.02338;
C_mu = 0;%0.05;?
C_malpha = -0.38;
C_malpha_dot = 0;%-6.7;?
C_mq = -3.6;
C_mTu =0;% -0.003;?
C_mTalpha = 0;%?

% Control Derivatives
C_Ddeltae = 0;
C_Dih = 0;%?
C_Ldeltae = -0.36;
C_Lih = 0;%0.94;%?
C_mdeltae = -0.5;
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
C_xu = -(C_Du+2*C_D1);
C_xalpha = -(C_Dalpha-C_L1);
C_xalpha_dot = 0;
C_xq = 0;
C_xdeltae = -C_Ddeltae;
C_xdeltap = C_Txu+2*C_Tx1;

C_zu = -(C_Lu+2*C_L1);
C_zalpha = -(C_Lalpha+C_D1);
C_zalpha_dot = -C_Lalpha_dot;
C_zq = -C_Lq;
C_zdeltae = -C_Ldeltae;
C_zdeltap = 0;

C_mu = C_mu+2*C_m1;
C_malpha = C_malpha;
C_malpha_dot = C_malpha_dot;
C_mq = C_mq;
C_mdeltae = C_mdeltae;
C_mdeltap = C_mTu+2*C_mT1;

% Nondimensional
C_w0 = w/(1/2*rho*u0^2*S);          %Non-dimensional weight [unitless]

% Longitudinal Dimensional Derivatives
Xu = rho*u0*S*C_w0*sin(theta0)+1/2*rho*u0*S*C_xu;
Xw = 1/2*rho*u0*S*C_xalpha;
Xq = 1/4*rho*u0*c_bar*S*C_xq;
Xw_dot = 1/4*rho*c_bar*S*C_xalpha_dot;

Zu = -rho*u0*S*C_w0*cos(theta0)+1/2*rho*u0*S*C_zu;
Zw = 1/2*rho*u0*S*C_zalpha;
Zq = 1/4*rho*u0*c_bar*S*C_zq;
Zw_dot = 1/4*rho*c_bar*S*C_zalpha_dot;

Mu = 1/2*rho*u0*c_bar*S*C_mu;
Mw = 1/2*rho*u0*c_bar*S*C_malpha;
Mq = 1/4*rho*u0*c_bar^2*S*C_mq;
Mw_dot = 1/4*rho*c_bar^2*S*C_malpha_dot;

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
X_deltae = C_xdeltae*1/2*rho*u0^2*S;
X_deltap = C_xdeltap*1/2*rho*u0^2*S;

Z_deltae = C_zdeltae*1/2*rho*u0^2*S;
Z_deltap = C_zdeltap*1/2*rho*u0^2*S;

M_deltae = C_mdeltae*1/2*rho*u0^2*S*c_bar;
M_deltap = C_mdeltap*1/2*rho*u0^2*S*c_bar;

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

Long1A = Xu/m;                                            %Row 1, Column 1
Long2A = Xw/m;                                            %Row 1, Column 2
Long3A = 0;                                               %Row 1, Column 3
Long4A = -g*cos(theta0);                                  %Row 1, Column 4
Long5A = Zu/(m-Zw_dot);                                   %Row 2, Column 1
Long6A = Zw/(m-Zw_dot);                                   %Row 2, Column 2
Long7A = (Zq + m*u0)/(m-Zw_dot);                          %Row 2, Column 3
Long8A = (-m*g*sin(theta0))/(m-Zw_dot);                   %Row 2, Column 4
Long9A = (1/Jy)*(Mu + (Mw_dot*Zu)/(m-Zw_dot));            %Row 3, Column 1
Long10A = (1/Jy)*(Mw + (Mw_dot*Zw)/(m-Zw_dot));           %Row 3, Column 2
Long11A = (1/Jy)*(Mq + (Mw_dot*(Zq+m*u0))/(m-Zw_dot));    %Row 3, Column 3
Long12A = -(Mw_dot*m*g*sin(theta0))/(Jy*(m-Zw_dot));      %Row 3, Column 4
Long13A = 0;                                              %Row 4, Column 1
Long14A = 0;                                              %Row 4, Column 2
Long15A = 1;                                              %Row 4, Column 3
Long16A = 0;                                              %Row 4, Column 4

A_long = [Long1A, Long2A, Long3A, Long4A;
          Long5A, Long6A, Long7A, Long8A;
          Long9A, Long10A, Long11A, Long12A;
          Long13A, Long14A, Long15A, Long16A]
      
B_long = [X_deltae/m                                     X_deltap/m;
          Z_deltae/(m-Zw_dot)                            Z_deltap/(m-Zw_dot);
          M_deltae/Jy+(Mw_dot*Z_deltae)/(Jy*(m-Zw_dot))  (M_deltap/Jy)+(Mw_dot*Z_deltap)/(Jy*(m-Zw_dot))
          0                                              0]
      
%% Lateral Linear Model
Lat1A = Yv/m;
Lat2A = Yp/m;
Lat3A = Yr/m - u0;
Lat4A = g*cos(theta0);
Lat5A = Lv/Ix_prime + Izx_prime*Nv;
Lat6A = Lp/Ix_prime + Izx_prime*Np;
Lat7A = Lr/Ix_prime + Izx_prime*Nr;
Lat8A = 0;
Lat9A = Izx_prime*Lv+Nv/Iz_prime;
Lat10A = Izx_prime*Lp+Np/Iz_prime;
Lat11A = Izx_prime*Lr+Nr/Iz_prime;
Lat12A = 0;
Lat13A = 0;
Lat14A = 1;
Lat15A = tan(theta0);
Lat16A = 0;

A_lat = [Lat1A, Lat2A, Lat3A, Lat4A;
          Lat5A, Lat6A, Lat7A, Lat8A;
          Lat9A, Lat10A, Lat11A, Lat12A;
          Lat13A, Lat14A, Lat15A, Lat16A]
      
B_lat = [Y_deltaa/m                              Y_deltar/m;
         (L_deltaa/Ix_prime)+Izx_prime*N_deltaa  (L_deltar/Ix_prime)+Izx_prime*N_deltar;
         (N_deltaa/Iz_prime)+Izx_prime*L_deltaa  (N_deltar/Iz_prime)+Izx_prime*L_deltar;
         0                                       0]

     
f_p_x = .5*rho*S_prop*C_prop*((k_motor*del_t)^2 - V_a^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs
% del_a
% del_p
% del_r_r
% del_r_l

[el_rud] = [1 1; -1 1]*[del_r_r; del_r_l];

del_e = el_rud(1);
del_r = el_rud(2);

X_long = [p_x;
          p_z;
          u_t;
          w_t;
          theta;
          q];

      
X_lat = [p_y;
         v_t;
         psi;
         phi;
         p;
         r];
     
U_long = [del_e;
          del_p];

U_lat = [del_a;
         del_r];
     
X_dot_long = A_long*X_long + B_long*U_long

X_dot_lat = A_lat*X_lat + B_lat*U_lat



f = [-m*g*sin(theta);
     m*g*cos(theta)*sin(phi);
     m*g*cos(theta)*cos(phi);] + 
     

