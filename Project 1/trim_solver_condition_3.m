function i = trim_solver_condition_3(f)
alpha_star = f(1);
beta_star = f(2);
phi_star = f(3);
% V_a_star : desired airspeed in[m/s]
% gamma_star : desired flight path angle [rad]
% R_star: desired turn radius
%% Aircraft Data (Aerosonde UAV)
% Geometric data
S = 0.55;                            %Wing surface area [m^2]
c = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]
%h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]
e = 0.9;                             %Oswald Efficiency Factor
%% Flight Conditions Data
g = 9.8;
altt = 100;                       %Altitude [m]
% Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
[T, a, P, rho_star] = atmosisa(altt); %Standard Atmospheric Calculations
% Mass and Inertial Data

m = 13.5;                            %Mass of the aircraft [kg]
W = m*g;                          %Weight of the aircraft [N]
J_x = 0.8244;                         %Moment of inertia x-axis [kg*m^2]
J_y = 1.135;                         %Moment of inertia y-axis [kg*m^2]
J_z = 1.759;                         %Moment of inertia z-axis [kg*m^2]
J_zx = 0.1204;                         %Moment of inertia zx-axis [kg*m^2]
J_xz = J_zx;                          %Moment of inertia xz-axis [kg*m^2]

Gamma = J_x*J_z - J_xz^2;

Gamma_1 = J_xz*(J_x-J_y + J_z)/Gamma;
Gamma_2 = (J_z*(J_z-J_y) + J_xz^2)/Gamma;
Gamma_3 = J_z/Gamma;
Gamma_4 = J_xz/Gamma;
Gamma_5 = (J_z - J_x)/J_y;
Gamma_6 = J_xz/J_y;
Gamma_7 = (J_x*(J_x-J_y) + J_xz^2)/Gamma;
Gamma_8 = J_x/Gamma;

%% Prop properties
S_prop = 0.2027; %m^2
k_motor = 80; 
C_prop = 1;

%% Known Trim Conditions
V_a_star = 30; %m/s
R_star = 150; %m
gamma_star = deg2rad(3); %rad

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
%% Aerodynamic Force Coefficients
%%% X direction
C_X = -C_D*cos(alpha_star)+ C_Lift*sin(alpha_star); %linear coefficient model
C_X_q = -C_D_q*cos(alpha_star)+ C_L_q*sin(alpha_star);
C_X_del_e = -C_D_delta_e*cos(alpha_star) + C_L_delta_e*sin(alpha_star);
%%% Z direction
C_Z = -C_D*sin(alpha_star) - C_Lift*cos(alpha_star);
C_Z_q = -C_D_q*sin(alpha_star)-C_L_q*cos(alpha_star);
C_Z_del_e = -C_D_delta_e*sin(alpha_star) - C_L_delta_e*cos(alpha_star);
%% Numerical Computation of Trim
%%%%% Solving for alpha_star, beta_star, phi_star
%%% Postion
% p_x_star = 0;
% p_y_star = 0;
% p_z_star = altt;
%%% Body frame velocities: u*,v*,w*
u_star = V_a_star*(cos(alpha_star)*cos(beta_star));
v_star = V_a_star*(sin(beta_star));
w_star = V_a_star*(sin(alpha_star)*cos(beta_star));
%%% Pitch angle (with Vw = 0)
theta_star = alpha_star + gamma_star;
%%% Angular rates (theta* is expressed in terms of gamma*and alpha*)
p_star = (V_a_star/R_star)*(-sin(theta_star));
q_star = (V_a_star/R_star)*(sin(phi_star)*cos(theta_star));
r_star = (V_a_star/R_star)*(cos(phi_star)*cos(theta_star));
%%% Elevator, del_e
delta_e_star = ((((J_xz*((p_star)^2 - (r_star)^2)+((J_x-J_z)*(p_star*r_star)))/(0.5*rho_star*(V_a_star)^2*c*S))) - C_m_0 - C_m_alpha*alpha_star - C_m_q*(c*q_star/(2*V_a_star)))/C_m_delta_e;
%%% Throttle, del_t
num1 = 2*m*(-r_star*v_star + q_star*w_star + g*sin(theta_star))-rho_star*(V_a_star)^2*S*(C_X + C_X_q*(c*q_star/(2*V_a_star))+ C_X_del_e*delta_e_star);
den1 = rho_star*S_prop*C_prop*(k_motor)^2;
delta_p_star = sqrt((num1/den1)+ ((V_a_star)^2/(k_motor)^2));
%%% Aileron (del_a) and Rudder (del_r)
m1 = [C_p_delta_a C_p_delta_r; C_r_delta_a C_r_delta_r];
m2_r1 = ((-Gamma_1*p_star*q_star + Gamma_2*q_star*r_star)/(0.5*rho_star*(V_a_star)^2*S*b))- C_p_0 - C_p_beta*beta_star - C_p_p*((b*p_star)/2*V_a_star) - C_p_r*((b*r_star/(2*V_a_star)));
m2_r2 = ((-Gamma_7*p_star*q_star + Gamma_1*q_star*r_star)/(0.5*rho_star*(V_a_star)^2*S*b))- C_r_0 - C_r_beta*beta_star - C_r_p*((b*p_star)/2*V_a_star) - C_r_r*((b*r_star/(2*V_a_star)));
m2 = [m2_r1;m2_r2];
a_r_star_matrix = inv(m1)*m2;
delta_a_star = a_r_star_matrix(1);
delta_r_star = a_r_star_matrix(2);
%%% f(x_star,u_star)
% F(1,1)= u_star;
% F(2,1) = v_star;
% F(3,1) = w_star;
% F(4,1) = theta_star;
% F(5,1) = p_star;
% F(6,1) = q_star;
% F(7,1) = r_star;
% F(8,1) = delta_e_star;
% F(9,1) = delta_t_star;
% F(10,1) = delta_a_star;
% F(11,1) = delta_r_star;
%% f(x*,u*)
%%% position_dot
% p_x_dot = (cos(theta_star)*cos(psi_star))*u_star+(sin(phi_star)*sin(theta)*cos(psi_star)- cos(phi_star)*sin(psi_star))*v_star +...
%     (cos(phi_star)*sin(theta_star)*cos(psi_star)+sin(phi_star)*sin(psi_star))*w_star;
% p_y_dot = (cos(theta_star)*sin(psi_star))*u_star+(sin(phi_star)*sin(theta_star)*sin(psi_star)+ cos(phi_star)*cos(psi_star))*v_star +...
%     (cos(phi_star)*sin(theta_star)*sin(psi_star)-sin(phi_star)*cos(psi_star))*w_star;
% p_z_dot = -1*(u_star*sin(theta_star)-v_star*sin(phi_star)*cos(theta_star)-w_star*cos(phi_star)*cos(theta_star)); %multiply by -1 b/c book defines h = -p_z
%%% wind_dot
u_dot = r_star*v_star - q_star*w_star - g*sin(theta_star)+ ...
    ((rho_star*(V_a_star)^2*S)/(2*m))*(C_X + C_X_q*((c*q_star)/(2*V_a_star))+ C_X_del_e*delta_e_star) +...
    ((rho_star*S_prop*C_prop)/(2*m))*((k_motor*delta_p_star)^2 - (V_a_star)^2);
v_dot = p_star*w_star - r_star*u_star + g*cos(theta_star)*sin(phi_star) + (rho_star*(V_a_star^2)*S/(2*m))*...
    (C_Y_0 + C_Y_beta*beta_star + C_Y_p*(b*p_star/(2*V_a_star))+ C_Y_r*(b*r_star/(2*V_a_star))+ C_Y_delta_a*delta_a_star + C_Y_delta_r*delta_r_star);
w_dot = q_star*u_star - p_star*v_star + g*cos(theta_star)*cos(phi_star) + (rho_star*(V_a_star^2)*S/(2*m))*(C_Z+ C_Z_q*(c*q_star/(2*V_a_star))+ C_Z_del_e*delta_e_star);
%%% euler rates
phi_dot = p_star + q_star*sin(phi_star)*tan(theta_star)+ r_star*cos(phi_star)*tan(theta_star);
theta_dot = q_star*cos(phi_star)- r_star*sin(phi_star);
% psi_dot = q_star*sin(phi_star)*sec(theta_star) + r_star*cos(phi_star)*sec(theta_star);
%%% angular rates
p_dot = Gamma_1*p_star*q_star - Gamma_2*q_star*r_star + 0.5*rho_star*(V_a_star^2)*S*b*...
    (C_p_0 + C_p_beta*beta_star + C_p_p*(b*p_star/(2*V_a_star))+ C_p_r*(b*r_star/(2*V_a_star))+ C_p_delta_a*delta_a_star + C_p_delta_r*delta_r_star);
q_dot = Gamma_5*p_star*r_star - Gamma_6*(p_star^2 - r_star^2) + (rho_star*(V_a_star^2)*S*c/(2*J_y))*...
    (C_m_0 + C_m_alpha*alpha_star + C_m_q*(c*q_star/(2*V_a_star))+ C_m_delta_e*delta_e_star);
r_dot = Gamma_7*p_star*q_star - Gamma_1*q_star*r_star + 0.5*rho_star*(V_a_star^2)*S*b*...
    (C_r_0 + C_r_beta*beta_star + C_r_p*(b*p_star/(2*V_a_star))+ C_r_r*(b*r_star/(2*V_a_star))+ C_r_delta_a*delta_a_star + C_r_delta_r*delta_r_star);

F(1,1) = u_dot;
F(2,1) = v_dot;
F(3,1) = w_dot;
F(4,1) = phi_dot;
F(5,1) = theta_dot;
% F(6,1) = psi_dot - (V_a_star/R_star)*cos(gamma_star);
F(6,1) = p_dot;
F(7,1) = q_dot;
F(8,1) = r_dot;

A = [0 0 0 0 0 0 0 0]';
i = norm(A-F)^2;
end