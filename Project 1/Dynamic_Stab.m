%% FROM ADC - B18 - HW3

clear all; clear all;
%% Constants needed for calculations (In Imperial units)
% g = 32.2;                           %Gravitational acceleration [ft/sec^2]
% R = 3089.2;                         %Specific gas constant for dry air [lb*ft/sl*K]
% a1 = -1.9812*10^-3;                 %Temperature gradient [K/ft]
% rho_sl = 0.0023769;                 %Density at sea level [sl/ft^3]
% tau_sl = 288.16;                    %Standard temperature at sea level [K]

%% Aircraft Data (Learjet 24 Aircraft)
% Geometric data
S = .55;                            %Wing surface area [m^2]
c_bar = 0.18994;                          %Mean Aerodynamic Chord (MAC) [m]
b = 2.8956;                             %Wing span [m]

% Flight Conditions Data
altt = 100;                       %Altitude [m]
%Mach = 0.7;                         %Mach of the aircraft [unitless]
V = 30;                            %True airspeed [m/sec]
h_CG = 0.32;                        %Esitmated location of the center of gravity on the wing chord [unitless]
alpha = deg2rad(2.7);               %Steady state angle of attck [rad]

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
C_Du = 0.104;?
C_Dalpha = 0.3;
C_Txu = -0.07;?
C_L0 = 0.28;
C_Lu = 0.4;?
C_Lalpha = 3.45;
C_Lalpha_dot = 2.2;?
C_Lq = 0;
C_m0 = -0.02338;
C_mu = 0.05;?
C_malpha = -0.38;
C_malpha_dot = -6.7;?
C_mq = -3.6;
C_mTu = -0.003;?
C_mTalpha = 0;?

% Control Derivatives
C_Ddeltae = 0;
C_Dih = 0;%?
C_Ldeltae = -0.36;
C_Lih = 0.94;%?
C_mdeltae = -0.5;
C_mih = -2.5;?

%% Lateral Directional Aerodynamic Coefficients
% Stability Derivatives
C_lbeta = -0.12;
C_lp = -0.26;
C_lr = 0.14;
C_Ybeta = -0.98;
C_Yp = 0;
C_Yr = 0;
C_nbeta = 0.25;
C_nTbeta = 0;?
C_np = 0.022;
C_nr = -0.35;

% Control Derivatives
C_ldeltaa = 0.08;
C_ldeltar = 0.105;
C_Ydeltaa = 0;
C_Ydeltar = -0.17;
C_ndeltaa = 0.06;
C_ndeltar = -0.074;?
%% Standard Atmospheric Calculations

[T, a, P, rho] = atmosisa(altt);

%     if altt <= 36089.2                                      %For altitudes 0 - 36089.2 km
%         tau = tau_sl + a1*altt;
%         rho = rho_sl * (tau/tau_sl).^(-1-g/(a1*R));
%     elseif altt > 36089.2                                   %For altitude 36089.2 - 82021 km
%         tau_11 = tau_sl + a1*36089.2;                       %Temperature at 36089.2 km
%         rho_11 = rho_sl * (tau_11/tau_sl).^(-1-g/(a1*R));   %Pressure at 36089.2 km
% 
%         tau = tau_11;
%         rho = rho_11*exp((-g*(altt-36089.2))/(R*tau_11));
%     end
    
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
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine if the aircraft is dynamically stable 
eig_long = eig(A_long)                 %Calculates the eigenvalues of A_long
eig_lat = eig(A_lat)                   %Calculates the eigenvalues of A_lat

% Longitudinal dynamic stability
if eig_long < 0
    fprintf('The aircraft is longitudinally dynamically stable \n')
else
    fprintf('The aircraft is not longitudinally dynamically stable \n')
end

% Lateral dynamic stability
if eig_lat < 0
    fprintf('The aircraft is laterally dynamically stable \n')
else
    fprintf('The aircraft is not laterally dynamically stable \n')
end

if eig_long < 0 & eig_lat < 0
    fprintf('The aircraft is dynamically stable \n')
else
    fprintf('The aircraft is not dynamically stable \n')
end
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine the settling times of the phugoid mode and of the roll convergence mode
%Phugoid mode
t_s_phugoid = 4/min(abs(real(eig_long)));
fprintf('The settling time for the phugoid mode is = %.3f seconds\n', t_s_phugoid)

%Roll convergence mode
t_s_roll_convg = 4/max(abs(real(eig_lat)));
fprintf('The settling time for the roll convergence mode is = %.3f seconds\n', t_s_roll_convg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the impulse and step response of the following input-output pairs: [elevator, pitch angle] [thrust, speed] [rudder, yaw] [aileron, roll angle]
    %All control inputs: [aileron (da), rudder (dr), elvator (de), thrust/throttle (dt)]
    %NOTE: longitudinal control inputs are: [elvator (de), thrust/throttle (dt)]
    %B_long matrix is set up such that "de" is the first input, and "dt" is the second input
    
long_sys = ss(A_long, B_long, eye(4), zeros(4,2));          %Creates a State System with state variables A, B, C, and D
figure; impulse(long_sys)                                   %Plots the impulse response for all input-output pairs
figure; step(long_sys)                                      %Plots the step response for all input-output pairs
lat_sys = ss([A_lat zeros(4,1); 0 0 sec(theta0) 0 0], [B_lat; 0 0], eye(5), zeros(5,2));
figure; impulse(lat_sys)
figure; step(lat_sys)

long_sys_el2pitch = ss(A_long, B_long(:,1)*pi/180, [0 0 0 1], 0);
long_sys_tr2speed = ss(A_long, B_long(:,2), [1 0 0 0], 0);
B_lat_wYaw = [B_lat; 0 0];

lat_sys_aileron2roll = ss(A_lat, B_lat(:,1)*pi/180, [0 0 0 1], 0);
%lat_sys_aileron2rollRate = ss(A_lat, B_lat(:,1)*pi/180, [0 1 0 0], 0);
lat_sys_rud2yaw = ss([A_lat zeros(4,1); 0 0 sec(theta0) 0 0], B_lat_wYaw(:,2)*pi/180, [0 0 0 0 1], 0);

figure
sgtitle("Elavator to Pitch")
subplot(121)
impulse(long_sys_el2pitch)
title("Impulse")
ylabel("\Delta\circ")
subplot(122)
step(long_sys_el2pitch)
title("Step")
ylabel("\Delta\circ")

figure
sgtitle("Aileron to roll")
subplot(121)
impulse(lat_sys_aileron2roll)
title("Impulse")
ylabel("\Delta\circ")
subplot(122)
step(lat_sys_aileron2roll)
title("Step")
ylabel("\Delta\circ")

figure
sgtitle("Thrust to Speed")
subplot(121)
impulse(long_sys_tr2speed)
title("Impulse")
ylabel("\Delta ft/s")
subplot(122)
step(long_sys_tr2speed)
title("Step")
ylabel("\Delta ft/s")

figure
sgtitle("Rudder to Yaw")
subplot(121)
impulse(lat_sys_rud2yaw)
title("Impulse")
ylabel("\Delta\circ")
subplot(122)
step(lat_sys_rud2yaw)
title("Step")
ylabel("\Delta\circ")

%% % Problem 2 % %% (Stick-fixed response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Controlled (closed-loop) response
% Create desired eigenvalues and thus settling times
desired_ev = [ -0.9873 + 2.6324i; -0.9873 - 2.6324i; -0.05 + 0.1115i; -0.05 - 0.1115i];
desired_settling_time = 4/0.05;

%Control design is about designing the right gain matrix K
K_long = place(A_long, B_long, desired_ev);

[t_sim_long, x_long] = ode45(@(t,x) (A_long-B_long*K_long)*x, [0 100], [10 0 0 0]');  %Last array is the distrubance introduced into the system

%% Rate of Change Plots %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta0 = 0;     %Initial pitch angle of the aircraft
u0 = 677;       %Initial velocity of the aircraft

A_long_with_altitude = [A_long, zeros(4,1);
                        -sin(theta0), cos(theta0), 0, -u0*cos(theta0), 0];
B_long_with_altitude = [B_long; 0 0];
                    
[t_sim_long, x_long] = ode45(@(t,xi) A_long_with_altitude*xi, [0 400], [10 0 0 0 0]');      %Stick-fixed response (no control component)

K_altitude_hold = lqr(A_long_with_altitude, B_long_with_altitude, 0.01*eye(5), 10*eye(2), zeros(5,2))

[t_sim_long_cl, x_long_cl] = ode45(@(t,xi_cl) (A_long_with_altitude-B_long_with_altitude*K_altitude_hold)*(xi_cl - [0; 0; 0; 0; 0]) , [0 400], [10 0 0 0 0]');

control_input_long = -K_altitude_hold*x_long_cl';

figure('Name', 'Longitudinal Stick-Fixed Response')
subplot(511)
plot(t_sim_long, x_long(:,1), t_sim_long_cl, x_long_cl(:,1))
xlabel('t (s)'); ylabel('\Delta u ft/s'); grid on;

subplot(512)
plot(t_sim_long, x_long(:,2), t_sim_long_cl, x_long_cl(:,2))
xlabel('t (s)'); ylabel('\Delta w ft/s'); grid on;

subplot(513)
plot(t_sim_long, x_long(:, 3)*180/pi, t_sim_long_cl, x_long_cl(:,3)*180/pi)
xlabel('t (s)'); ylabel('\Delta q (\circ/s)'); grid on;

subplot(514)
plot(t_sim_long, x_long(:, 4)*180/pi, t_sim_long_cl, x_long_cl(:,4)*180/pi)
xlabel('t (s)'); ylabel('\Delta \theta (\circ)'); grid on;

subplot(515)
plot(t_sim_long, x_long(:,5), t_sim_long_cl, x_long_cl(:,5))
xlabel('t (s)'); ylabel('\Delta z_E ft'); grid on; 

figure
plot(t_sim_long_cl, control_input_long(1,:)*180/pi')
title("Elevator control input over time")
xlabel("t(s)"); ylabel("\delta \circ")
xlim([0 5])
figure
plot(t_sim_long_cl, control_input_long(2,:)*180/pi')
title("Thrust control input over time")
xlabel("t(s)"); ylabel("\Delta lbs")

%% Lateral stick-fixed response
%%% Part 2 %%% - Roll Angle
[t_sim_lat, x_sim_lat] = ode45(@(t,x) A_lat*x, [0 100], [0 1*pi/180 0 0]');

K_lat = lqr(A_lat, B_lat, 0.01*eye(4), 10*eye(2), zeros(4,2));

[t_sim_cl, x_sim_lat_cl] = ode45(@(t,x) (A_lat - B_lat*K_lat)*(x - [0; 0; 0; 5*pi/180]) , [0 100], [0 1*pi/180 0 0]');

control_input_lat = -K_lat*x_sim_lat_cl';

figure('Name', 'Lateral stick fixed');
subplot(411)
plot(t_sim_lat, x_sim_lat(:, 1), t_sim_cl, x_sim_lat_cl(:,1))
xlabel('t (s)'); ylabel('v ft/s'); grid on;

subplot(412)
plot(t_sim_lat, x_sim_lat(:, 2)*180/pi,  t_sim_cl, x_sim_lat_cl(:,2)*180/pi)
xlabel('t (s)'); ylabel('p (\circ/s)'); grid on;


subplot(413)
plot(t_sim_lat, x_sim_lat(:, 3)*180/pi,  t_sim_cl, x_sim_lat_cl(:,3)*180/pi)
xlabel('t (s)'); ylabel('r (\circ/s)'); grid on;


subplot(414)
plot(t_sim_lat, x_sim_lat(:, 4)*180/pi,  t_sim_cl, x_sim_lat_cl(:,4)*180/pi)
xlabel('t (s)'); ylabel('\phi (\circ)'); grid on; 

figure
plot(t_sim_cl, control_input_lat*180/pi')
title("Lateral control Inputs over time")
xlabel("t(s)"); ylabel("(\circ)")
legend('\delta a','\delta r')
xlim([0 10]);

figure
sideSlip = asin(x_sim_lat_cl(:,1)/u0);
plot(t_sim_cl, sideSlip*180/pi);
title("Side Slip angle over time")
xlabel("t(s)"); ylabel(" \beta (\circ)")
xlim([0 5]);