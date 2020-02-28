clear all; clc;
%% Aircraft Data (Aerosonde UAV)
% Geometric data
S = 0.55;                            %Wing surface area [m^2]
c = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]
%h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]
e = 0.9;
%% Flight Conditions Data

g = 9.8;
altt = 100;                       %Altitude [m]
Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
[T, a, P, rho_star] = atmosisa(altt); %Standard Atmospheric Calculations
%alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

%% Trim

V_a_star = 30; %m/s
R_star = inf;
gamma_star = 0;
theta_star = deg2rad(4); 
psi_star = 0;
phi_star = 0; 
u_star = 0;
v_star = 0;
w_star = 0;
p_star = 0;
q_star = 0;
r_star = 0;
% Mass and Inertial Data

m = 13.5;                            %Mass of the aircraft [kg]
W = m*g;                          %Weight of the aircraft [N]
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
C_L_0 = 0.28;
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
%% Control Inputs and Angles
%%% AOA and Sideslip
alpha_and_delta_e = inv([C_L_alpha C_L_delta_e; C_m_alpha C_m_delta_e])*...
    [(W/(0.5*rho_star*V^2*S)-C_L_0); -C_m_0];
alpha_star = alpha_and_delta_e(1);
alpha_ig = alpha_star;
beta_star = asin(v_star/V_a_star);
%%% Control Inputs
delta_e_star = alpha_and_delta_e(2);
delta_a_star = 0; 
delta_r_star = 0;
% delta_p* calculation
AR = (b^2)/S; %aspect ratio
k = (4/3)*(1/(pi*e*AR)); %characterizes induced drag
C_Lift = W/(0.5*rho_star*(V^2)*S);
L_D_ratio = C_Lift/(C_D_0 + k*(C_Lift)^2);
T = W/L_D_ratio; %thrust
delta_p_star = (1/k_motor)*sqrt(((2*T)/(rho_star*S_prop*C_prop))+ (V_a_star^2));

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

C_X_0 = -C_D_0*cos(alpha_star) + C_L_0*sin(alpha_star);
C_X_alpha = -C_D_alpha*cos(alpha_star) + C_L_alpha*sin(alpha_star);
C_X_q = -C_D_q*cos(alpha_star) + C_L_q*sin(alpha_star);
C_X_delta_e = -C_D_delta_e*cos(alpha_star) + C_L_delta_e*sin(alpha_star);

C_Z_0 = -C_D_0*sin(alpha_star) + C_L_0*cos(alpha_star);
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
%% Coefficient of Lift and Drag
C_Lift = C_L_0 + C_L_alpha*alpha_star;
C_D = C_D_0 + C_D_alpha*alpha_star;
L = 0.5*rho_star*(V_a_star^2)*S*C_Lift;
%% Aerodynamic Force Coefficients
%%% X direction
C_X = -C_D*cos(alpha_star)+ C_Lift*sin(alpha_star); %linear coefficient model
C_X_q = -C_D_q*cos(alpha_star)+ C_L_q*sin(alpha_star);
C_X_del_e = -C_D_delta_e*cos(alpha_star) + C_L_delta_e*sin(alpha_star);
%%% Z direction
C_Z = -C_D*sin(alpha_star) - C_Lift*cos(alpha_star);
C_Z_q = -C_D_q*sin(alpha_star)-C_L_q*cos(alpha_star);
C_Z_del_e = -C_D_delta_e*sin(alpha_star) - C_L_delta_e*cos(alpha_star);

%% Longitudinal Dimensional Derivatives
Xu = ((u_star*rho_star*S)/m)*((C_X_0 + C_X_alpha*alpha_star + C_X_delta_e*delta_e_star) - ((rho_star*S*w_star*C_X_alpha)/(2*m)) + ((rho_star*S*c*C_X_q*u_star*q_star)/4*m*V_a_star) - ((rho_star*S_prop*C_prop*u_star)/m));
Xw = -q_star + ((w_star*rho_star*S)/m) * (C_X_0 + C_X_alpha*alpha_star + C_X_delta_e*delta_e_star) + ((rho_star*S*c*C_X_q*w_star*q_star)/4*m*V_a_star) + ((rho_star*S*u_star*C_X_alpha)/2*m) - ((rho_star*S_prop*C_prop*w_star)/m);
Xq = -w_star + ((rho_star*V_a_star*c*S*C_X_q)/4*m);

Zu = q_star + ((u_star*rho_star*S)/m) * (C_Z_0 + C_Z_alpha*alpha_star + C_Z_delta_e*delta_e_star) - ((rho_star*S*w_star*C_Z_alpha)/2*m) + ((rho_star*S*c*C_Z_q*u_star*q_star)/4*m*V_a_star);
Zw = ((w_star*rho_star*S)/m) * (C_Z_0 + C_Z_alpha*alpha_star + C_Z_delta_e*delta_e_star) - ((rho_star*S*u_star*C_Z_alpha)/2*m) + ((rho_star*S*c*C_Z_q*w_star*q_star)/4*m*V_a_star);
Zq = u_star + ((rho_star*V_a_star*S*C_Z_q*q_star)/4*m);

Mu = ((u_star*rho_star*S*c)/J_y)*(C_m_0 + C_m_alpha*alpha_star + C_m_delta_e*delta_e_star) - ((rho_star*S*c*C_m_alpha*u_star)/2*J_y) + (rho_star*S*c^2 *q_star*u_star)/(4*J_y*V_a_star);
Mw = ((w_star*rho_star*S*c)/J_y)*(C_m_0 + C_m_alpha*alpha_star + C_m_delta_e*delta_e_star) + ((rho_star*S*c*C_m_alpha*u_star)/2*J_y) + (rho_star*S*c^2 *q_star*u_star)/(4*J_y*V_a_star);
Mq = ((1/4*rho_star*V_a_star*c^2*S*C_m_q)/J_y);

%Lateral Dimensional Derivatives
Yv = ((rho_star*S*b*v_star)/(4*m*V_a_star))*(C_Y_p*p_star + C_Y_r*r_star) + ((rho_star*S*v_star)/m)*(C_Y_0 + C_Y_beta*beta_star + C_Y_delta_a*delta_a_star + C_Y_delta_r*delta_r_star) + ((rho_star*S*C_Y_beta)/(2*m))*sqrt(u_star^2+w_star^2);
Yp = w_star + (rho_star*V_a_star*S*b)/(4*m)*C_Y_p;
Yr = u_star + (rho_star*V_a_star*S*b)/(4*m)*C_Y_r;

Lv = ((rho_star*S*b*v_star)/(4*V_a_star))*(C_p_p*p_star + C_p_r*r_star) + (rho_star*S*v_star)*(C_p_0 + C_p_beta*beta_star + C_p_delta_a*delta_a_star + C_p_delta_r*delta_r_star) + ((rho_star*S*C_p_beta)/(2))*sqrt(u_star^2+w_star^2);
Lp = Gamma_1*q_star+ (rho_star*V_a_star*S*b^2)/4 * C_p_p;
Lr = -Gamma_2*q_star+ (rho_star*V_a_star*S*b^2)/4 * C_p_r;

Nv = ((rho_star*S*b*v_star)/(4*V_a_star))*(C_r_p*p_star + C_r_r*r_star) + (rho_star*S*v_star)*(C_r_0 + C_r_beta*beta_star + C_r_delta_a*delta_a_star + C_r_delta_r*delta_r_star) + ((rho_star*S*C_r_beta)/2)*sqrt(u_star^2+w_star^2);
Np = Gamma_7*q_star+ (rho_star*V_a_star*S*b^2)/4 * C_r_p;
Nr = -Gamma_1*q_star+ (rho_star*V_a_star*S*b^2)/4 * C_r_r;

%% Dimensional Control Derivatives
%Longitudinal Dimensional Control Derivatives
X_delta_e = (C_X_delta_e*1/2*rho_star*V_a_star^2*S)/m;
X_delta_t = (rho_star*S_prop*C_prop*k_motor^2*delta_p_star)/m;

Z_delta_e = (C_Z_delta_e*1/2*rho_star*V_a_star^2*S)/m;

M_delta_e = (C_m_delta_e*1/2*rho_star*V_a_star^2*S*c)/J_y;

%Lateral Dimensional Control Derivatives
Y_delta_a = (C_Y_delta_a*1/2*rho_star*V_a_star^2*S)/m;
Y_delta_r = (C_Y_delta_r*1/2*rho_star*V_a_star^2*S)/m;

L_delta_a = C_p_delta_a*1/2*rho_star*V_a_star^2*S*b;
L_delta_r = C_p_delta_r*1/2*rho_star*V_a_star^2*S*b;

N_delta_a = C_r_delta_a*1/2*rho_star*V_a_star^2*S*b;
N_delta_r = C_r_delta_r*1/2*rho_star*V_a_star^2*S*b;
%% Finding Trim States
%% trim_solver_condition_1
% "Constant speed, constant altitude cruise"
V_a_star = 30; %m/s
R_star = inf; %m
gamma_star = 0; %rad
%fsolve
IGs = [alpha_ig 0 0];
trim_alpha_beta_phi = fsolve(@trim_solver_condition_1, IGs);
alpha_star = trim_alpha_beta_phi(1)
beta_star = trim_alpha_beta_phi(2)
phi_star = trim_alpha_beta_phi(3);
% trim states
ts = trim_states(V_a_star,R_star,gamma_star,alpha_star,beta_star,phi_star);
u_star = ts(1);
v_star = ts(2);
w_star = ts(3);
theta_star = ts(4);
p_star = ts(5);
q_star = ts(6);
r_star = ts(7);
trim_states_x_star_1 = [u_star, v_star, w_star, theta_star, p_star, q_star, r_star]';
%trim control inputs
tci = trim_control_inputs(alpha_star,beta_star,g,V_a_star,rho_star,J_x,J_z,J_xz,c,m,b,S,k_motor,S_prop,C_prop,v_star,w_star,C_m_0, C_m_alpha,C_m_q,C_m_delta_e,C_p_r,C_r_r,C_X,C_X_q,C_X_del_e,theta_star,p_star,q_star,r_star,C_p_delta_a,C_p_delta_r, C_r_delta_a,C_r_delta_r,C_p_0,C_r_0,C_p_p,C_r_p,C_p_beta,C_r_beta,Gamma_1,Gamma_2,Gamma_7);
delta_e_star = tci(1);
delta_p_star = tci(2);
delta_a_star = tci(3);
delta_r_star = tci(4);
trim_control_inputs_1 = [delta_e_star, delta_p_star, delta_a_star, delta_r_star]';
fprintf('Trim Condition: 1')
trim_state_x_star = [u_star v_star w_star theta_star p_star q_star r_star]'
trim_controlinputs = [delta_e_star delta_p_star delta_a_star delta_r_star]'
[A_long,B_long] = longitudinal_linear_model(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, X_delta_e, X_delta_t,Z_delta_e, M_delta_e, g, u_star, w_star, theta_star);
[A_lat,B_lat] = lateral_linear_model(Yv,Yp,Yr,Lv,Lp,Lr,Nv,Np,Nr,Y_delta_a,Y_delta_r,L_delta_a,L_delta_r,N_delta_a,N_delta_r,theta_star,phi_star,p_star,q_star,r_star,g);
Along_eigenvectors = eig(A_long)
Alat_eigenvectors = eig(A_lat)
%% trim_solver_condition_2
% "Constant speed, constant altitude turn"
V_a_star = 30; %m/s
R_star = 150; %m
gamma_star = 0; %rad
%fsolve
IGs = [alpha_ig 0 0];
trim_alpha_beta_phi = fminunc(@trim_solver_condition_2, IGs);
alpha_star = trim_alpha_beta_phi(1);
beta_star = trim_alpha_beta_phi(2);
phi_star = trim_alpha_beta_phi(3);
% trim states
ts = trim_states(V_a_star,R_star,gamma_star,alpha_star,beta_star,phi_star);
u_star = ts(1);
v_star = ts(2);
w_star = ts(3);
theta_star = ts(4);
p_star = ts(5);
q_star = ts(6);
r_star = ts(7);
trim_states_x_star_2 = [u_star, v_star, w_star, theta_star, p_star, q_star, r_star]';
%trim control inputs
tci = trim_control_inputs(alpha_star,beta_star,g,V_a_star,rho_star,J_x,J_z,J_xz,c,m,b,S,k_motor,S_prop,C_prop,v_star,w_star,C_m_0, C_m_alpha,C_m_q,C_m_delta_e,C_p_r,C_r_r,C_X,C_X_q,C_X_del_e,theta_star,p_star,q_star,r_star,C_p_delta_a,C_p_delta_r, C_r_delta_a,C_r_delta_r,C_p_0,C_r_0,C_p_p,C_r_p,C_p_beta,C_r_beta,Gamma_1,Gamma_2,Gamma_7);
delta_e_star = tci(1);
delta_p_star = tci(2);
delta_a_star = tci(3);
delta_r_star = tci(4);
trim_control_inputs_2 = [delta_e_star, delta_p_star, delta_a_star, delta_r_star]';
fprintf('Trim Condition: 2')
trim_state_x_star = [u_star v_star w_star theta_star p_star q_star r_star]'
trim_controlinputs = [delta_e_star delta_p_star delta_a_star delta_r_star]'
[A_long,B_long] = longitudinal_linear_model(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, X_delta_e, X_delta_t,Z_delta_e, M_delta_e, g, u_star, w_star, theta_star);
[A_lat,B_lat] = lateral_linear_model(Yv,Yp,Yr,Lv,Lp,Lr,Nv,Np,Nr,Y_delta_a,Y_delta_r,L_delta_a,L_delta_r,N_delta_a,N_delta_r,theta_star,phi_star,p_star,q_star,r_star,g);
Along_eigenvectors = eig(A_long)
Alat_eigenvectors = eig(A_lat)
%% trim_solver_condition_3
% "Constant speed climbing turn"
V_a_star = 30; %m/s
R_star = 150; %m
gamma_star = deg2rad(3); %rad
%fsolve
% phi_ig = asin(R_star*L/(m*V_a_star^2));
IGs = [alpha_ig 0 0];
trim_alpha_beta_phi = fminunc(@trim_solver_condition_3, IGs);
alpha_star = trim_alpha_beta_phi(1);
beta_star = trim_alpha_beta_phi(2);
phi_star = trim_alpha_beta_phi(3);
% trim states
ts = trim_states(V_a_star,R_star,gamma_star,alpha_star,beta_star,phi_star);
u_star = ts(1);
v_star = ts(2);
w_star = ts(3);
theta_star = ts(4);
p_star = ts(5);
q_star = ts(6);
r_star = ts(7);
trim_states_x_star_3 = [u_star, v_star, w_star, theta_star, p_star, q_star, r_star]';
%trim control inputs
tci = trim_control_inputs(alpha_star,beta_star,g,V_a_star,rho_star,J_x,J_z,J_xz,c,m,b,S,k_motor,S_prop,C_prop,v_star,w_star,C_m_0, C_m_alpha,C_m_q,C_m_delta_e,C_p_r,C_r_r,C_X,C_X_q,C_X_del_e,theta_star,p_star,q_star,r_star,C_p_delta_a,C_p_delta_r, C_r_delta_a,C_r_delta_r,C_p_0,C_r_0,C_p_p,C_r_p,C_p_beta,C_r_beta,Gamma_1,Gamma_2,Gamma_7);
delta_e_star = tci(1);
delta_p_star = tci(2);
delta_a_star = tci(3);
delta_r_star = tci(4);
trim_control_inputs_3 = [delta_e_star, delta_p_star, delta_a_star, delta_r_star]';fprintf('Trim Condition: 3')
trim_state_x_star = [u_star v_star w_star theta_star p_star q_star r_star]'
trim_controlinputs = [delta_e_star delta_p_star delta_a_star delta_r_star]'
[A_long,B_long] = longitudinal_linear_model(Xu, Xw, Xq, Zu, Zw, Zq, Mu, Mw, Mq, X_delta_e, X_delta_t,Z_delta_e, M_delta_e, g, u_star, w_star, theta_star);
[A_lat,B_lat] = lateral_linear_model(Yv,Yp,Yr,Lv,Lp,Lr,Nv,Np,Nr,Y_delta_a,Y_delta_r,L_delta_a,L_delta_r,N_delta_a,N_delta_r,theta_star,phi_star,p_star,q_star,r_star,g);
Along_eigenvectors = eig(A_long)
Alat_eigenvectors = eig(A_lat)



