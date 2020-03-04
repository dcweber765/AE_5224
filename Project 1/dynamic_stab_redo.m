clear all

S = .55;                            %Wing surface area [m^2]
c_bar = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]
%h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]

%% Flight Conditions Data

g = 9.8;
altt = 100;                       %Altitude [m]
Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
[T, a, P, rho_star] = atmosisa(altt); %Standard Atmospheric Calculations
%alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

% AR = b^2/S;
% lamda = 0;
% LAMDA = 0;
% dEpsilon = 4.44*(sqrt(1-Mach^2))*(((1/AR) - (1/(1+AR^1.7)))*((10-3*lamda)/7)*(1-(ltv/b)/(2*(lt/b)^.33))*sqrt(cos(LAMDA)))^1.19;
%% Trim

V_star = 30;
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
theta_star = 0;                %Initial pitch angle [rad]
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
C_D_u = 0;
C_L_0 = 0.28;
C_L_alpha = 3.45;
C_L_q = 0;
C_L_u = 0;
C_M_0 = -0.02338;
C_M_alpha = -0.38;
C_M_q = -3.6;
C_M_u = 0;

C_D_p = 0.0437;

C_D_star = C_D_0+(C_D_alpha*alpha_star);
C_L_star = C_L_0+(C_L_alpha*alpha_star);
C_M_alpha_dot = C_M_q/3; %Pg 259 [6]
C_L_alpha_dot = 0; %guess

% Control Derivatives
C_D_delta_e = 0;
C_L_delta_e = -0.36;
C_M_delta_e = -0.5;

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


Q_star_S = .5*rho_star*V_star*S;

X_u = 2*((m*g)/V_star)*sin(theta_star) - Q_star_S*(C_D_u + 2*C_D_star);
X_w = Q_star_S*(C_L_star - C_D_alpha);
Z_u = -2*((m*g)/V_star)*cos(theta_star) - Q_star_S*(C_L_u + 2*C_L_star);
Z_w = -Q_star_S*(C_D_star + C_L_alpha);
Z_q = -.5*Q_star_S*c_bar*C_L_q;
Z_w_dot = (-Q_star_S/(2*V_star))*c_bar*C_L_alpha_dot;
M_u = Q_star_S*c_bar*C_M_u;
M_w = Q_star_S*c_bar*C_M_alpha;
M_q = .5*Q_star_S*c_bar^2*C_M_q;
M_w_dot = (Q_star_S/(2*V_star))*c_bar^2*C_M_alpha_dot;

X_delta_e = -Q_star_S*V_star*C_D_delta_e;
Z_delta_e = -Q_star_S*V_star*C_L_delta_e;
M_delta_e = Q_star_S*V_star*c_bar*C_M_delta_e;

X_delta_t = rho_star*S_prop*C_prop*k_motor^2*delta_t_star;
Z_delta_t = 0;
%M_delta_t = Q_star_S*V_star*c_bar*C_M_delta_t;


M_prime = M_w_dot/(m - Z_w_dot);

ALong_1_1 = X_u/m;
ALong_1_2 = X_w/m;
ALong_1_3 = 0;
ALong_1_4 = -g*cos(theta_star);
ALong_1_5 = 0;
ALong_1_6 = 0;

ALong_2_1 = Z_u/(m - Z_w_dot);
ALong_2_2 = Z_w/(m - Z_w_dot);
ALong_2_3 = (Z_q + m*u_star)/(m - Z_w_dot);
ALong_2_4 = -m*g*sin(theta_star)/(m - Z_w_dot);
ALong_2_5 = 0;
ALong_2_6 = 0;

ALong_3_1 = (M_u + M_prime*Z_u)/J_y;
ALong_3_2 = (M_w + M_prime*Z_w)/J_y;
ALong_3_3 = (M_q + M_prime*(Z_q+m*u_star))/J_y;
ALong_3_4 = (-M_prime*m*g*sin(theta_star))/J_y;
ALong_3_5 = 0;
ALong_3_6 = 0;

ALong_4_1 = 0;
ALong_4_2 = 0;
ALong_4_3 = 1;
ALong_4_4 = 0;
ALong_4_5 = 0;
ALong_4_6 = 0;

ALong_5_1 = -sin(theta_star);
ALong_5_2 = cos(theta_star);
ALong_5_3 = 0;
ALong_5_4 = -V_star*cos(theta_star);
ALong_5_5 = 0;
ALong_5_6 = 0;

ALong_6_1 = cos(theta_star);
ALong_6_2 = sin(theta_star);
ALong_6_3 = 0;
ALong_6_4 = -V_star*sin(theta_star);
ALong_6_5 = 0;
ALong_6_6 = 0;

BLong_1_1 = X_delta_e/m;
BLong_1_2 = X_delta_t/m;

BLong_2_1 = Z_delta_e/(m - Z_w_dot);
BLong_2_2 = Z_delta_t/(m - Z_w_dot);

BLong_3_1 = M_delta_e/J_y + ((M_w_dot*Z_delta_e)/(J_y*(m - Z_w_dot)));
BLong_3_2 = 0;%C_M_delta_t/J_y + ((M_w_dot*Z_delta_t)/(J_y*(m - Z_w_dot)));

BLong_4_1 = 0;
BLong_4_2 = 0;

BLong_5_1 = 0;
BLong_5_2 = 0;
 
BLong_6_1 = 0;
BLong_6_2 = 0;


ALong = [ALong_1_1 ALong_1_2 ALong_1_3 ALong_1_4 ALong_1_5 ALong_1_6;...
         ALong_2_1 ALong_2_2 ALong_2_3 ALong_2_4 ALong_2_5 ALong_2_6;...
         ALong_3_1 ALong_3_2 ALong_3_3 ALong_3_4 ALong_3_5 ALong_3_6;...
         ALong_4_1 ALong_4_2 ALong_4_3 ALong_4_4 ALong_4_5 ALong_4_6;...
         ALong_5_1 ALong_5_2 ALong_5_3 ALong_5_4 ALong_5_5 ALong_5_6;...
         ALong_6_1 ALong_6_2 ALong_6_3 ALong_6_4 ALong_6_5 ALong_6_6]

BLong = [BLong_1_1 BLong_1_2;...
         BLong_2_1 BLong_2_2;...  
         BLong_3_1 BLong_3_2;...  
         BLong_4_1 BLong_4_2;...  
         BLong_5_1 BLong_5_2;...  
         BLong_6_1 BLong_6_2];  
     
Long_sys = ss(ALong,BLong,eye(6),0);     
x0 = [.1;0;0;0;10;0]
%t = linspace(0,300,600);
%u  = [delta_e_star*ones(1,600); delta_t_star*ones(1,600)]; %delta_e_star delta_t_star
figure; initial(Long_sys,x0,200)
%lsim(Long_sys, u, t, x0)
