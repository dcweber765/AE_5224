function x_lat_dot = lat_dynamics(x, un, star)

delta_a = un(1);
delta_r = un(2);

if delta_a > deg2rad(20)
    delta_a = deg2rad(20);
elseif delta_a < deg2rad(-20)
    delta_a = deg2rad(-20);
end


if delta_r > deg2rad(20)
    delta_r = deg2rad(20);
elseif delta_r < deg2rad(-20)
    delta_r = deg2rad(-20);
end

un(1) = delta_a;
un(2) = delta_r;

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
C_N_0 = 0;
C_N_beta = 0.25;
C_N_p = 0.022;
C_N_r = -0.35;

% Control Derivatives
C_l_delta_a = 0.08;
C_l_delta_r = 0.105;
C_Y_delta_a = 0;
C_Y_delta_r = -0.17;
C_N_delta_a = 0.06;
C_N_delta_r = -0.032;


Q_star_S = .5*rho_star*V_star*S;

Y_v = Q_star_S*C_Y_beta;
Y_p = .5*Q_star_S*C_Y_p;
Y_r = .5*Q_star_S*b*C_Y_r;
l_v = Q_star_S*b*C_l_beta;
l_p = Q_star_S*b^2*C_l_p;
l_r = .5*Q_star_S*b^2*C_l_r;
N_v = Q_star_S*b*C_N_beta;
N_p = .5*Q_star_S*b^2*C_N_p;
N_r = .5*Q_star_S*b^2*C_N_r;

Y_delta_a = Q_star_S*V_star*C_Y_delta_a;
Y_delta_r = Q_star_S*V_star*C_Y_delta_r;

l_delta_a = Q_star_S*V_star*b*C_l_delta_a;
l_delta_r = Q_star_S*V_star*b*C_l_delta_r;

N_delta_a = Q_star_S*V_star*b*C_N_delta_a;
N_delta_r = Q_star_S*V_star*b*C_N_delta_r;

J_x_prime = (J_x*J_z-J_zx^2)/J_z;
J_z_prime = (J_x*J_z-J_zx^2)/J_x;
J_zx_prime = J_zx/(J_x*J_z-J_xz^2);

ALat_1_1 = Y_v/m;
ALat_1_2 = Y_p/m;
ALat_1_3 = -V_star + Y_r/m;
ALat_1_4 = g*cos(theta_star);
ALat_1_5 = 0;
ALat_1_6 = 0;

ALat_2_1 = l_v/J_x_prime + J_zx_prime*N_v;
ALat_2_2 = l_p/J_x_prime + J_zx_prime*N_p;
ALat_2_3 = l_r/J_x_prime + J_zx_prime*N_r;
ALat_2_4 = 0;
ALat_2_5 = 0;
ALat_2_6 = 0;

ALat_3_1 = N_v/J_z_prime + J_zx_prime*l_v;
ALat_3_2 = N_p/J_z_prime + J_zx_prime*l_p;
ALat_3_3 = N_r/J_z_prime + J_zx_prime*l_r;
ALat_3_4 = 0;
ALat_3_5 = 0;
ALat_3_6 = 0;

ALat_4_1 = 0;
ALat_4_2 = 1;
ALat_4_3 = tan(theta_star);
ALat_4_4 = 0;
ALat_4_5 = 0;
ALat_4_6 = 0;

ALat_5_1 = 0;
ALat_5_2 = 0;
ALat_5_3 = sec(theta_star);
ALat_5_4 = 0;
ALat_5_5 = 0;
ALat_5_6 = 0;

ALat_6_1 = 1;
ALat_6_2 = 0;
ALat_6_3 = 0;
ALat_6_4 = 0;
ALat_6_5 = V_star*cos(theta_star);
ALat_6_6 = 0;


BLat_1_1 = Y_delta_a/m;
BLat_1_2 = Y_delta_r/m;

BLat_2_1 = l_delta_a/J_x_prime + J_zx_prime*N_delta_a;
BLat_2_2 = l_delta_r/J_x_prime + J_zx_prime*N_delta_r;

BLat_3_1 = N_delta_a/J_z_prime + J_zx_prime*l_delta_a;
BLat_3_2 = N_delta_r/J_z_prime + J_zx_prime*l_delta_r;

BLat_4_1 = 0;
BLat_4_2 = 0;

BLat_5_1 = 0;
BLat_5_2 = 0;

BLat_6_1 = 0;
BLat_6_2 = 0;

ALat = [ALat_1_1 ALat_1_2 ALat_1_3 ALat_1_4 ALat_1_5 ALat_1_6;...
    ALat_2_1 ALat_2_2 ALat_2_3 ALat_2_4 ALat_2_5 ALat_2_6;...
    ALat_3_1 ALat_3_2 ALat_3_3 ALat_3_4 ALat_3_5 ALat_3_6;...
    ALat_4_1 ALat_4_2 ALat_4_3 ALat_4_4 ALat_4_5 ALat_4_6;...
    ALat_5_1 ALat_5_2 ALat_5_3 ALat_5_4 ALat_5_5 ALat_5_6;...
    ALat_6_1 ALat_6_2 ALat_6_3 ALat_6_4 ALat_6_5 ALat_6_6];

BLat = [BLat_1_1 BLat_1_2;...
    BLat_2_1 BLat_2_2;...
    BLat_3_1 BLat_3_2;...
    BLat_4_1 BLat_4_2;...
    BLat_5_1 BLat_5_2;...
    BLat_6_1 BLat_6_2];

x_lat_dot = ALat*x + BLat*un;

end
