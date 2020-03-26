function x_long_dot = long_dynamics(x, un, star)

V_star = star(1);
alpha_star = star(2);
beta_star = star(3);
delta_e_star = star(4);
delta_t_star = star(5);
delta_a_star = star(6);
delta_r_star = star(7);

u_star = star(8);
v_star = star(9);
w_star = star(10);

p_star = star(11);
q_star = star(12);
r_star = star(13);

psi_star = star(14);
theta_star = star(15);
phi_star = star(16);

m = star(17);

J_x = star(18);
J_y = star(19);
J_z = star(20);
J_zx = star(21);
J_xz = star(22); 

rho_star = star(23);
S = star(24);
c_bar = star(25);
b =star(26);
g = 9.81;

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

%% Long
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
BLong_3_2 = 0;%M_delta_t/J_y + ((M_w_dot*Z_delta_t)/(J_y*(m - Z_w_dot)));

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
    ALong_6_1 ALong_6_2 ALong_6_3 ALong_6_4 ALong_6_5 ALong_6_6];

BLong = [BLong_1_1 BLong_1_2;...
    BLong_2_1 BLong_2_2;...
    BLong_3_1 BLong_3_2;...
    BLong_4_1 BLong_4_2;...
    BLong_5_1 BLong_5_2;...
    BLong_6_1 BLong_6_2];


x_long_dot = ALong*x + BLong*un;
end