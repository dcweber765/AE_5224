clear variables
close all

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
    ALong_6_1 ALong_6_2 ALong_6_3 ALong_6_4 ALong_6_5 ALong_6_6]

BLong = [BLong_1_1 BLong_1_2;...
    BLong_2_1 BLong_2_2;...
    BLong_3_1 BLong_3_2;...
    BLong_4_1 BLong_4_2;...
    BLong_5_1 BLong_5_2;...
    BLong_6_1 BLong_6_2];

%% Lat

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
    ALat_6_1 ALat_6_2 ALat_6_3 ALat_6_4 ALat_6_5 ALat_6_6]

BLat = [BLat_1_1 BLat_1_2;...
    BLat_2_1 BLat_2_2;...
    BLat_3_1 BLat_3_2;...
    BLat_4_1 BLat_4_2;...
    BLat_5_1 BLat_5_2;...
    BLat_6_1 BLat_6_2];

eigLong = eig(ALong)
eigLat = eig(ALat)

Long_sys = ss(ALong,BLong,eye(6),0);
Lat_sys = ss(ALat,BLat,eye(6),0);
x0_Long = [.01;0;0;0;.1;0];
x0_Lat = [0;0;0;0;.01;0];
%t = linspace(0,300,600);
%u  = [delta_e_star*ones(1,600); delta_t_star*ones(1,600)]; %delta_e_star delta_t_star
figure; initial(Long_sys,x0_Long)
figure; initial(Lat_sys,x0_Lat)
%lsim(Long_sys, u, t, x0)
[t_sim_long, x_long] = ode45(@(t,x) ALong*x, [0 400], [1 0 0 0 0 0]');
[t_sim_lat, x_lat] = ode45(@(t,x) ALat*x, [0 400], [0 1*pi/180 0 0 0 0]');


figure('Name', 'Longitudinal Stick-Fixed Response')
subplot(511)
plot(t_sim_long, x_long(:,1))
xlabel('t (s)'); ylabel('\Delta u ft/s'); grid on;

subplot(512)
plot(t_sim_long, x_long(:,2))
xlabel('t (s)'); ylabel('\Delta w ft/s'); grid on;

subplot(513)
plot(t_sim_long, x_long(:, 3)*180/pi)
xlabel('t (s)'); ylabel('\Delta q (\circ/s)'); grid on;

subplot(514)
plot(t_sim_long, x_long(:, 4)*180/pi)
xlabel('t (s)'); ylabel('\Delta \theta (\circ)'); grid on;

subplot(515)
plot(t_sim_long, x_long(:,5))
xlabel('t (s)'); ylabel('\Delta z_E ft'); grid on;


figure('Name', 'Lateral stick fixed');
subplot(511)
plot(t_sim_lat, x_lat(:, 1))
xlabel('t (s)'); ylabel('v ft/s'); grid on;

subplot(512)
plot(t_sim_lat, x_lat(:, 2)*180/pi)
xlabel('t (s)'); ylabel('p (\circ/s)'); grid on;


subplot(513)
plot(t_sim_lat, x_lat(:, 3)*180/pi)
xlabel('t (s)'); ylabel('r (\circ/s)'); grid on;


subplot(514)
plot(t_sim_lat, x_lat(:, 4)*180/pi)
xlabel('t (s)'); ylabel('\phi (\circ)'); grid on;

subplot(515)
plot(t_sim_lat, x_lat(:,5))
xlabel('t (s)'); ylabel('\Delta z_E ft'); grid on;



% p_x = x_long(:,5);
% p_y = x_lat(:,6);
% p_z = x_long(:,6);
% u_b = x_long(:,1);
% v_b = x_lat(:,1);
% w_b = x_long(:,2);
% u_t = 0;
% v_t = 0;
% w_t = 0;
% V = sqrt(u_t.^2 .* v_t.^2 .* w_t^2);
% alpha = atan(w_b./u_b);
% beta = asin(v_b./V);
% psi = x_lat(:,5);
% theta = x_long(:,4);
% phi = x_lat(:,4);
% p = x_lat(:,2);
% q = x_long(:,3);
% r = x_lat(:,3);
% 
% % p_t = [p_x p_y p_z];
% % V_t = [u_t v_t w_t];
% % V_b = [u_b v_b w_b];
% % euler_ang = [psi theta phi];
% % omega_tb_b = [p q r];

%X = [p_t V_t V_b V euler_ang omega_tb_b]' %alpha;beta;
n_pts	= 1001;
time_pts= linspace(0, 300, 1001);
xLong_from_true	= zeros(6, n_pts);
control_in_long = zeros(2,n_pts);
x_long_true =	x0_Long;
xLong_from_true(:, 1)= x_long_true;
starT1 = [V_star,alpha_star,beta_star,delta_e_star,delta_t_star,delta_a_star,delta_r_star,u_star,v_star,w_star,p_star,q_star,r_star,psi_star,theta_star,phi_star,m,J_x,J_y,J_z,J_zx,J_xz, rho_star,S,c_bar,b];
%% LONG
for m1 = 2:n_pts
    
    % RK4 step (Euler)
    dt	= time_pts(m1) - time_pts(m1 - 1);
    u1	= control_in_long(:, m1 - 1);
    u2	= control_in_long(:, m1);
    u12	= 0.5*(u1 + u2);							% Approx (p,q,r) at (t + 0.5dt)
    
    k1	= dt*long_dynamics(x_long_true,				u1, starT1);
    k2	= dt*long_dynamics((x_long_true + 0.5*k1),	u12, starT1);
    k3	= dt*long_dynamics((x_long_true + 0.5*k2),	u12, starT1);
    k4	= dt*long_dynamics((x_long_true + k3),		u2, starT1);
    x_long_true	= x_long_true + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
    
    % Record
    xLong_from_true(:, m1)	= x_long_true;
end
%% LAT
xLat_from_true	= zeros(6, n_pts);
control_in_lat = zeros(2,n_pts);
x_lat_true =	x0_Lat;
xLat_from_true(:, 1)= x_lat_true;
for m1 = 2:n_pts
    
    % RK4 step (Euler)
    dt	= time_pts(m1) - time_pts(m1 - 1);
    u1	= control_in_lat(:, m1 - 1);
    u2	= control_in_lat(:, m1);
    u12	= 0.5*(u1 + u2);							% Approx (p,q,r) at (t + 0.5dt)
    
    k1	= dt*lat_dynamics(x_lat_true,				u1, starT1);
    k2	= dt*lat_dynamics((x_lat_true + 0.5*k1),	u12, starT1);
    k3	= dt*lat_dynamics((x_lat_true + 0.5*k2),	u12, starT1);
    k4	= dt*lat_dynamics((x_lat_true + k3),		u2, starT1);
    x_lat_true	= x_lat_true + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
    
    % Record
    xLat_from_true(:, m1)	= x_lat_true;
end

p_x = xLong_from_true(5,:);
p_y = xLat_from_true(6,:);
p_z = xLong_from_true(6,:);
u_b = xLong_from_true(1,:);
v_b = xLat_from_true(1,:);
w_b = xLong_from_true(2,:);
u_t = xLong_from_true(1,:);
v_t = xLat_from_true(1,:);
w_t = xLong_from_true(2,:);
V = sqrt(u_t.^2 .* v_t.^2 .* w_t.^2);
alpha = atan(w_b./u_b);
beta = asin(v_b./V);
psi = xLat_from_true(5,:);
theta = xLong_from_true(4,:);
phi = xLat_from_true(4,:);
p = xLat_from_true(2,:);
q = xLong_from_true(3,:);
r = xLat_from_true(3,:);

p_t = [p_x; p_y; p_z];
V_t = [u_t; v_t; w_t];
V_b = [u_b; v_b; w_b];
euler_ang = [psi; theta; phi];
omega_tb_b = [p; q; r];

X = [p_t; V_t; V_b; V; euler_ang; omega_tb_b;];

figure('Name', 'Longitudinal Stick-Fixed Response RK4')
subplot(511)
plot(time_pts, xLong_from_true(1,:))
xlabel('t (s)'); ylabel('\Delta u ft/s'); grid on;

subplot(512)
plot(time_pts, xLong_from_true(2,:))
xlabel('t (s)'); ylabel('\Delta w ft/s'); grid on;

subplot(513)
plot(time_pts, xLong_from_true(3,:)*180/pi)
xlabel('t (s)'); ylabel('\Delta q (\circ/s)'); grid on;

subplot(514)
plot(time_pts, xLong_from_true(3,:)*180/pi)
xlabel('t (s)'); ylabel('\Delta \theta (\circ)'); grid on;

subplot(515)
plot(time_pts, xLong_from_true(5,:))
xlabel('t (s)'); ylabel('\Delta z_E ft'); grid on;


figure('Name', 'Lateral stick fixed RK4');
subplot(511)
plot(time_pts, xLat_from_true(1,:))
xlabel('t (s)'); ylabel('v ft/s'); grid on;

subplot(512)
plot(time_pts, xLat_from_true(2,:)*180/pi)
xlabel('t (s)'); ylabel('p (\circ/s)'); grid on;


subplot(513)
plot(time_pts, xLat_from_true(3,:)*180/pi)
xlabel('t (s)'); ylabel('r (\circ/s)'); grid on;


subplot(514)
plot(time_pts, xLat_from_true(4,:)*180/pi)
xlabel('t (s)'); ylabel('\phi (\circ)'); grid on;

subplot(515)
plot(time_pts, xLat_from_true(5,:))
xlabel('t (s)'); ylabel('\Delta z_E ft'); grid on;
