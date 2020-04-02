close all; clear variables;

%% Aircraft geometry

S = .55;                            %Wing surface area [m^2]
c_bar = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]

%% Flight Conditions Data

g = 9.8;
altt = 100;                       %Altitude [m]
Mach = 0.1;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
[T, a, P, rho_star] = atmosisa(altt); %Standard Atmospheric Calculations

%% Trim conditons

V_star = 30;
alpha_star = 0;
beta_star = -.2737;
delta_e_star = -.0448;
delta_t_star = .4190;
delta_a_star = .7242;
delta_r_star = -.8764;

u_star = 28.8836;                             %Initial velocity [m/s]
v_star = -8.1078;
w_star = 0;

p_star = 0;
q_star = -.0817;
r_star = .1825;

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


x0_Long = [.01;0;0;0;.1;0];
x0_Lat = [0;0;0;0;.01;0];

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
subplot(611)
plot(time_pts, xLong_from_true(1,:))
xlabel('t (s)'); ylabel('\Delta u m/s'); grid on;

subplot(612)
plot(time_pts, xLong_from_true(2,:))
xlabel('t (s)'); ylabel('\Delta w m/s'); grid on;

subplot(613)
plot(time_pts, xLong_from_true(3,:)*180/pi)
xlabel('t (s)'); ylabel('\Delta q (\circ/s)'); grid on;

subplot(614)
plot(time_pts, xLong_from_true(4,:)*180/pi)
xlabel('t (s)'); ylabel('\Delta \theta (\circ)'); grid on;

subplot(615)
plot(time_pts, xLong_from_true(5,:))
xlabel('t (s)'); ylabel('\Delta p_z m'); grid on;

subplot(616)
plot(time_pts, xLong_from_true(6,:))
xlabel('t (s)'); ylabel('\Delta p_x m'); grid on;


figure('Name', 'Lateral stick fixed RK4');
subplot(611)
plot(time_pts, xLat_from_true(1,:))
xlabel('t (s)'); ylabel('v ft/s'); grid on;

subplot(612)
plot(time_pts, xLat_from_true(2,:)*180/pi)
xlabel('t (s)'); ylabel('p (\circ/s)'); grid on;


subplot(613)
plot(time_pts, xLat_from_true(3,:)*180/pi)
xlabel('t (s)'); ylabel('r (\circ/s)'); grid on;


subplot(614)
plot(time_pts, xLat_from_true(4,:)*180/pi)
xlabel('t (s)'); ylabel('\phi (\circ)'); grid on;

subplot(615)
plot(time_pts, xLat_from_true(5,:)*180/pi)
xlabel('t (s)'); ylabel('\psi (\circ)'); grid on;

subplot(616)
plot(time_pts, xLat_from_true(6,:))
xlabel('t (s)'); ylabel('\Delta p_y ft'); grid on;