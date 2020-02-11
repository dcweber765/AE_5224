function localization_ekf()
close all; clc;

%% Synthetic measurements
n_pts	= 1001;
time_pts= linspace(0, 300, 1001);	% s

generating_polynomial = 0.5*[
    0.3192   -0.0301    1.0933    0.0774   -0.0068    0.3714   -1.0891;
    0.3129   -0.1649    1.1093   -1.2141    1.5326   -0.2256    0.0326;
   -0.8649    0.6277   -0.8637   -1.1135   -0.7697    1.1174    0.5525];

euler_angles_true		= zeros(3, n_pts);
quaternions_true		= zeros(4, n_pts);
euler_angles_true_dot	= zeros(3, n_pts);
pqr_true				= zeros(3, n_pts);
for m1 = 1:3
	euler_angles_true(m1, :)= ( generating_polynomial(m1, 1)*(sin(0.1*time_pts)) + ...
		generating_polynomial(m1, 2)*(sin(0.2*time_pts)) + ...
		generating_polynomial(m1, 3)*(sin(0.3*time_pts)) + ...
		generating_polynomial(m1, 4)*(cos(0.03*time_pts)) + ...
		generating_polynomial(m1, 5)*(cos(0.02*time_pts)) + ...
		generating_polynomial(m1, 6)*(cos(0.5*time_pts)) + generating_polynomial(m1, 7) )/2;
	
	euler_angles_true_dot(m1, :) = ( 0.1*generating_polynomial(m1, 1)*(cos(0.1*time_pts)) + ...
		0.2*generating_polynomial(m1, 2)*(cos(0.2*time_pts)) + ...
		0.3*generating_polynomial(m1, 3)*(cos(0.3*time_pts)) - ...
		0.03*generating_polynomial(m1, 4)*(sin(0.03*time_pts)) - ...
		0.02*generating_polynomial(m1, 5)*(sin(0.02*time_pts)) - ...
		0.5*generating_polynomial(m1, 6)*(sin(0.5*time_pts)) )/2;
end


for m1 = 1:n_pts
	e_psi	= euler_angles_true(1, m1);
	e_thta	= euler_angles_true(2, m1);
	e_phi	= euler_angles_true(3, m1);
	pqr_true(:, m1) = [-sin(e_thta) 0 1; ...
			sin(e_phi)*cos(e_thta) cos(e_phi) 0; ...
			cos(e_phi)*cos(e_thta) -sin(e_phi) 0] * euler_angles_true_dot(:, m1);
		
	quaternions_true(:, m1)	= euler2quat( euler_angles_true(:, m1) );
	
end

std_dev_mag = 0.3*pi/180;					% Magnetometer std. dev. (rad)
std_dev_gyro= 0.13*pi/180;					% Rate gyro std. dev. (rad/s)
std_dev_airspeed= 0.002;                    % Airspeed sensor std. dev (kPa)
std_dev_accel= 0.0025*9.81;                 % Accelerometer std. dev (m/s^2)
std_dev_gps= 0.21%*9.009e-6;                 % GPS std. dev (deg)
std_dev_gps_alt= 0.4%*9.009e-6;              % GPS std. dev for altitude (deg)
std_dev_gps_vel= 0.05;                      % GPS velocity std. dev (m/s)
bias_mag	= 1*pi/180 * randn(3, 1);       % Magnetometer bias (rad)
bias_gyro	= 1*pi/180 * randn(3, 1);       % Rate gyro bias (rad/s)
bias_airspeed= 0.020;                       % Airspeed sensor bias (kPa)
bias_gps= 4.7%*9.009e-6;                     % GPS sensor bias (deg)
bias_accel= 0.01*9.81;                      % Accelerometer sensor bias (m/s^2) *guessed*
z_euler_angles_noisy	= euler_angles_true;
z_rate_gyro_noisy		= pqr_true;
for m1 = 1:numel(time_pts)
	z_euler_angles_noisy(:, m1)	= euler_angles_true(:, m1) + std_dev_mag*randn(3, 1) + bias_mag;
	z_rate_gyro_noisy(:, m1)	= pqr_true(:, m1) + std_dev_gyro*randn(3, 1) + bias_gyro;
end



save attitude_kinematics_synthetic_data.mat time_pts pqr_true ...
	euler_angles_true z_euler_angles_noisy z_rate_gyro_noisy ...
	std_dev_gyro std_dev_mag

% subplot(311); plot(time_pts, euler_angles_true(1,:)*180/pi, time_pts, z_euler_angles_noisy(1,:)*180/pi); xlabel('Time (s)'); ylabel('\psi (deg)'); grid on;
% subplot(312); plot(time_pts, euler_angles_true(2,:)*180/pi, time_pts, z_euler_angles_noisy(2,:)*180/pi); xlabel('Time (s)'); ylabel('\theta (deg)');  grid on;
% subplot(313); plot(time_pts, euler_angles_true(3,:)*180/pi, time_pts, z_euler_angles_noisy(3,:)*180/pi); xlabel('Time (s)'); ylabel('\phi (deg)');  grid on;

%% "Streaming" integration of Euler angles and quaternions (true)

euler_angles_int_from_true	= zeros(3, n_pts);
quaternions_int_from_true	= zeros(4, n_pts);
quaternions_int_from_true_norm	= zeros(1, n_pts);

x_e =	euler_angles_true(:, 1);
x_q	=	quaternions_true(:, 1);
euler_angles_int_from_true(:, 1)= x_e;
quaternions_int_from_true(:, 1)	= x_q;
for m1 = 2:n_pts
	
	% RK4 step (Euler)
	dt	= time_pts(m1) - time_pts(m1 - 1);
	u1	= pqr_true(:, m1 - 1);
	u2	= pqr_true(:, m1);
	u12	= 0.5*(u1 + u2);							% Approx (p,q,r) at (t + 0.5dt)
	
	k1	= dt*euler_kinematics(x_e,				u1);
	k2	= dt*euler_kinematics((x_e + 0.5*k1),	u12);
	k3	= dt*euler_kinematics((x_e + 0.5*k2),	u12);
	k4	= dt*euler_kinematics((x_e + k3),		u2);
	x_e	= x_e + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
	
	k1	= dt*quaternion_kinematics(x_q,				u1);
	k2	= dt*quaternion_kinematics((x_q + 0.5*k1),	u12);
	k3	= dt*quaternion_kinematics((x_q + 0.5*k2),	u12);
	k4	= dt*quaternion_kinematics((x_q + k3),		u2);
	x_q	= x_q + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
	
	% Record
	euler_angles_int_from_true(:, m1)	= x_e;
	quaternions_int_from_true(:, m1)	= x_q;
	quaternions_int_from_true_norm(:, m1)= norm(x_q);
end

%% "Streaming" integration of Euler angles and quaternions (noisy, unfiltered)

euler_angles_int_from_noisy = zeros(3, n_pts);
quaternions_int_from_noisy	= zeros(4, n_pts);
quaternions_int_from_noisy_norm	= zeros(1, n_pts);

x_e =	z_euler_angles_noisy(:, 1);
x_q =	euler2quat(x_e);

euler_angles_int_from_noisy(:, 1)	= x_e;
quaternions_int_from_noisy(:, 1)	= x_q;
for m1 = 2:n_pts
	
	% RK4 step
	dt	= time_pts(m1) - time_pts(m1 - 1);
	u1	= z_rate_gyro_noisy(:, m1 - 1);
	u2	= z_rate_gyro_noisy(:, m1);
	u12	= 0.5*(u1 + u2);							% Approx (p,q,r) at (t + 0.5dt)
	k1	= dt*euler_kinematics(x_e,				u1);
	k2	= dt*euler_kinematics((x_e + 0.5*k1),	u12);
	k3	= dt*euler_kinematics((x_e + 0.5*k2),	u12);
	k4	= dt*euler_kinematics((x_e + k3),		u2);
	x_e	= x_e + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
	
	k1	= dt*quaternion_kinematics(x_q,				u1);
	k2	= dt*quaternion_kinematics((x_q + 0.5*k1),	u12);
	k3	= dt*quaternion_kinematics((x_q + 0.5*k2),	u12);
	k4	= dt*quaternion_kinematics((x_q + k3),		u2);
	x_q	= x_q + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;

% 	x_q = x_q + dt*quaternion_kinematics(x_q, u1)
	
	% Record
	euler_angles_int_from_noisy(:, m1) = x_e;
	quaternions_int_from_noisy(:, m1)	= x_q;
	quaternions_int_from_noisy_norm(:, m1)= norm(x_q);
end

%% "Streaming" extended Kalman filter estimates of Euler angles and biases

z_gps_pos = ones(3,n_pts);
z_gps_velocity = ones(3,n_pts);
z_accel_noisy = ones(3,n_pts);

R_me_covariance		= blkdiag((std_dev_gps^2)*eye(2), (std_dev_gps_alt^2),(std_dev_gps_vel^2)*eye(3),(std_dev_mag^2) * eye(3));								% Measurement error covariance (3x3 matrix)
P_ee_covariance     = diag(1000*ones(18,1));						

% Initial estimate same as direct measurement; zero initial estimate of all biases
x_e_hat	= [z_gps_pos(:,1); z_gps_velocity(:,1); z_euler_angles_noisy(:,1); zeros(9,1)]; %bias_mag; bias_mag; bias_mag; bias_accel; bias_accel; bias_accel; bias_gyro; bias_gyro; bias_gyro];
x_e_dot = [z_gps_velocity(:,1); z_accel_noisy(:,1); z_rate_gyro_noisy(:,1); zeros(9,1)];


euler_angles_ekf_est	= zeros(3, n_pts);
% biases_e_ekf_est		= zeros(9, n_pts);		% (three magnetometer, three rate gyro) biases
% trace_P_rec				= zeros(1, n_pts);

euler_angles_ekf_est(:, 1)	= x_e_hat(7:9);
biases_e_ekf_est(:, 1)		= zeros(6, 1);		% (three magnetometer, three rate gyro) biases
trace_P_rec(1)				= trace(P_ee_covariance); 

for m1 = 2:n_pts
    dt = time_pts(m1) - time_pts(m1 - 1);
	
	% Physics-based predictive model 
    a_x = z_accel_noisy(1);
    a_y = z_accel_noisy(2);
    a_z = z_accel_noisy(3);
    z_p = z_rate_gyro_noisy(1);
    z_q = z_rate_gyro_noisy(2);
    z_r = z_rate_gyro_noisy(3);
    e_psi   = x_e_hat(7);
	e_thta	= x_e_hat(8);
	e_phi	= x_e_hat(9);
    Hinv	= (1 / cos(e_thta) ) * [0 sin(e_phi) cos(e_phi); ...
		0 cos(e_phi)*cos(e_thta) -sin(e_phi)*cos(e_thta); ...
		cos(e_thta) sin(e_phi)*sin(e_thta) cos(e_phi)*sin(e_thta)];
    
    R_t_b = [cos(e_psi)*cos(e_thta) cos(e_thta)*sin(e_psi) -sin(e_thta);...
        cos(e_psi)*sin(e_phi)*sin(e_thta)-cos(e_phi)*sin(e_psi) cos(e_phi)*cos(e_psi)+sin(e_phi)*sin(e_psi)*sin(e_thta) cos(e_thta)*sin(e_phi);...
        sin(e_phi)*sin(e_psi)+cos(e_phi)*cos(e_psi)*sin(e_thta) cos(e_phi)*sin(e_psi)*sin(e_thta)-cos(e_psi)*sin(e_phi) cos(e_phi)*cos(e_thta)];
    
    Jac_Hinv = [0 z_q*sin(e_phi)*sin(e_thta)/(cos(e_thta)^2)+z_r*cos(e_phi)*sin(e_thta)/(cos(e_thta)^2) z_q*cos(e_phi)/cos(e_thta)-z_r*sin(e_phi)/cos(e_thta);...
        0 0 -z_q*sin(e_phi)-z_r*cos(e_phi);...
        0 z_q*sin(e_phi)/(cos(e_thta)^2)+z_r*cos(e_phi)/(cos(e_thta)^2) z_q*cos(e_phi)*tan(e_thta)-z_r*sin(e_phi)*tan(e_thta)];
    
    Jac_Rtb = [-a_x*sin(e_psi)*cos(e_thta)+a_y*cos(e_thta)*cos(e_psi) -a_x*cos(e_psi)*sin(e_thta)-a_y*sin(e_thta)*sin(e_psi)-a_z*cos(e_thta) 0;...
        a_x*(-sin(e_psi)*sin(e_phi)*sin(e_thta)-cos(e_phi)*cos(e_psi))+a_y*(sin(e_phi)*cos(e_psi)*sin(e_thta)-cos(e_phi)*sin(e_psi)) a_x*cos(e_psi)*sin(e_phi)*cos(e_thta)+a_y*sin(e_phi)*sin(e_psi)*cos(e_thta)-a_z*sin(e_phi)*sin(e_thta) a_x*(cos(e_psi)*cos(e_phi)*sin(e_thta)+sin(e_phi)*sin(e_psi))+a_y*(cos(e_phi)*sin(e_psi)*sin(e_thta)-sin(e_phi)*cos(e_psi))+a_z*cos(e_thta)*cos(e_phi);...
        a_x*(sin(e_phi)*cos(e_psi)-cos(e_phi)*sin(e_psi)*sin(e_thta))+a_y*(cos(e_phi)*cos(e_psi)*sin(e_thta)+sin(e_psi)*sin(e_phi)) a_x*cos(e_phi)*cos(e_psi)*cos(e_thta)+a_y*cos(e_phi)*sin(e_psi)*cos(e_thta)-a_z*cos(e_phi)*sin(e_thta) a_x*(cos(e_phi)*sin(e_psi)-sin(e_phi)*cos(e_psi)*sin(e_thta))+a_y*(-sin(e_phi)*sin(e_psi)*sin(e_thta)-cos(e_psi)*cos(e_phi))-a_z*sin(e_phi)*cos(e_thta)];
	
	A_e	= [zeros(3) eye(3) zeros(3,12);...
        zeros(3) zeros(3) Jac_Rtb zeros(3) -R_t_b zeros(3);...
        zeros(3) zeros(3) Jac_Hinv zeros(3,6) -Hinv;...
        zeros(9,18)];
	B_e	= [zeros(3,6);...
        -R_t_b zeros(3);...
        zeros(3) -Hinv;...
        zeros(9,6)];
	Q_pn_covariance = B_e * (blkdiag((std_dev_accel^2)*eye(3),(std_dev_gyro^2) * eye(3))) * B_e';
	
	zm	= [z_gps_pos(:,m1); z_gps_velocity(:,m1); z_euler_angles_noisy(:,m1); zeros(9,1)];	% New magnetometer measurements
	
	% Predictive estimate (RK4) and e.e. covariance
% 	u1	= z_rate_gyro_noisy(:, m1 - 1) - biases_e_ekf_est(4:6, m1 - 1);
% 	u2	= z_rate_gyro_noisy(:, m1) - x_e_hat(7:9);
% 	u12	= 0.5*(u1 + u2);													% Approx (p,q,r) at (t + 0.5dt)
% 	k1	= dt*euler_kinematics(x_e_hat(1:3),				u1);
% 	k2	= dt*euler_kinematics((x_e_hat(1:3) + 0.5*k1),	u12);
% 	k3	= dt*euler_kinematics((x_e_hat(1:3) + 0.5*k2),	u12);
% 	k4	= dt*euler_kinematics((x_e_hat(1:3) + k3),		u2);
%  	x_e_minus	= x_e_hat + [ ...
%  		((1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4); zeros(6,1)];			% The predictive state update can use nonlinear kinematic equations

    x_e_minus	= x_e_hat + x_e_dot*dt;
    
    P_minus= (eye(18) + A_e*dt)*P_ee_covariance*(eye(18) + A_e*dt)' + Q_pn_covariance*dt^2;			% The e.e. covariance update uses linearized A...
	
	% Measurement model
    C_e	= eye(18);		% For Euler angles
	
	% Kalman gain
    K_e_gain = P_minus*C_e' * inv(C_e*P_minus*C_e' + blkdiag(R_me_covariance,0*eye(9)))			%... as does the Kalman gain computation
	
	% Recursive WLSE method to update predictive estimate with new
	% measurement
	x_e_hat	= x_e_minus + K_e_gain * ( zm - C_e*x_e_minus );				% New estimate
	P_ee_covariance		= (eye(18) - K_e_gain * C_e) * P_minus	;		% New e.e. covariance, again linearized C 
	
	% Record
% 	euler_angles_ekf_est(:, m1)			= x_e_hat(7:9);
% 	biases_e_ekf_est(:, m1)				= x_e_hat(10:end);
% 	trace_P_rec(m1)						= trace(P_ee_covariance);
	% Notice that the estimation error covariance decreases with time, 
	% indicating increasing confidence in the position estimate
end


figsize = [0, 0.04, 0.8, 0.7];
figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);

subplot(321); plot(time_pts, pqr_true(1,:)*180/pi, 'LineWidth', 2); hold on;
subplot(321); plot(time_pts, z_rate_gyro_noisy(1,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('\omega_{tb}^{b} (deg/s)'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(323); plot(time_pts, pqr_true(2,:)*180/pi, 'LineWidth', 2); hold on;
subplot(323); plot(time_pts, z_rate_gyro_noisy(2,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('\omega_{tb}^{b} (deg/s)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(325); plot(time_pts, pqr_true(3,:)*180/pi, 'LineWidth', 2); hold on;
subplot(325); plot(time_pts, z_rate_gyro_noisy(3,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('\omega_{tb}^{b} (deg/s)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(322); plot(time_pts, euler_angles_true(1,:)*180/pi, 'LineWidth', 2); hold on;
subplot(322); plot(time_pts, z_euler_angles_noisy(1,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('z_{\psi} (deg)'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(324); plot(time_pts, euler_angles_true(2,:)*180/pi, 'LineWidth', 2); hold on;
subplot(324); plot(time_pts, z_euler_angles_noisy(2,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('z_{\theta} (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(326); plot(time_pts, euler_angles_true(3,:)*180/pi, 'LineWidth', 2); hold on;
subplot(326); plot(time_pts, z_euler_angles_noisy(3,:)*180/pi, 'LineWidth', 1); hold on;
xlabel('Time (s)'); ylabel('z_{\phi} (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

print -dpng attitude_ekf_data_constbias.png



figsize = [0, 0.04, 0.5, 0.5];
figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);
plot(time_pts, trace_P_rec, 'LineWidth', 2); grid on;
xlabel('Time (s)'); ylabel('tr(P)');
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
print -dpng attitude_ekf_traceP_constbias.png


figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);

subplot(311); plot(time_pts, euler_angles_int_from_true(1,:)*180/pi, 'LineWidth', 2); hold on;
subplot(311); plot(time_pts, euler_angles_int_from_noisy(1,:)*180/pi, 'LineWidth', 1); 
subplot(311); plot(time_pts, euler_angles_ekf_est(1,:)*180/pi, '--', 'LineWidth', 2); 
xlabel('Time (s)'); ylabel('\psi (deg)'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(312); plot(time_pts, euler_angles_int_from_true(2,:)*180/pi, 'LineWidth', 2); hold on;
subplot(312); plot(time_pts, euler_angles_int_from_noisy(2,:)*180/pi, 'LineWidth', 1);
subplot(312); plot(time_pts, euler_angles_ekf_est(2,:)*180/pi, '--', 'LineWidth', 2); 
xlabel('Time (s)'); ylabel('\theta (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(313); plot(time_pts, euler_angles_int_from_true(3,:)*180/pi, 'LineWidth', 2); hold on;
subplot(313); plot(time_pts, euler_angles_int_from_noisy(3,:)*180/pi, 'LineWidth', 1);
subplot(313); plot(time_pts, euler_angles_ekf_est(3,:)*180/pi, '--', 'LineWidth', 2); 
xlabel('Time (s)'); ylabel('\phi (deg)');  grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

print -dpng attitude_ekf_euler_constbias.png

figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);
for m1 = 1:6
	subplot(2, 3, m1);
	plot(time_pts, biases_e_ekf_est(m1, :)*180/pi, 'LineWidth', 2);
	set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
	hold on; grid on;
end
print -dpng attitude_ekf_est_bias.png

figsize = [0, 0.04, 0.5, 0.7];
figure('Units', 'Normalized', 'InnerPosition', figsize, 'OuterPosition', figsize);
% subplot(511); plot(time_pts, quaternions_true(1,:), 'LineWidth', 2); hold on;
subplot(511); plot(time_pts, quaternions_int_from_true(1,:), 'LineWidth', 2); hold on;
subplot(511); plot(time_pts, quaternions_int_from_noisy(1,:), 'LineWidth', 1);
xlabel('Time (s)'); ylabel('e_0'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

% subplot(512); plot(time_pts, quaternions_true(2,:), 'LineWidth', 2); hold on;
subplot(512); plot(time_pts, quaternions_int_from_true(2,:), 'LineWidth', 2); hold on;
subplot(512); plot(time_pts, quaternions_int_from_noisy(2,:), 'LineWidth', 1);
xlabel('Time (s)'); ylabel('e_1'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

% subplot(513); plot(time_pts, quaternions_true(3,:), 'LineWidth', 2); hold on;
subplot(513); plot(time_pts, quaternions_int_from_true(3,:), 'LineWidth', 2); hold on;
subplot(513); plot(time_pts, quaternions_int_from_noisy(3,:), 'LineWidth', 1);
xlabel('Time (s)'); ylabel('e_2'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

% subplot(514); plot(time_pts, quaternions_true(4,:), 'LineWidth', 2); hold on;
subplot(514); plot(time_pts, quaternions_int_from_true(4,:), 'LineWidth', 2); hold on;
subplot(514); plot(time_pts, quaternions_int_from_noisy(4,:), 'LineWidth',1);
xlabel('Time (s)'); ylabel('e_3'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');

subplot(515); plot(time_pts, quaternions_int_from_true_norm, 'LineWidth', 2); hold on;
subplot(515); plot(time_pts, quaternions_int_from_noisy_norm, 'LineWidth', 1);
xlabel('Time (s)'); ylabel('||e||'); grid on;
set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
print -dpng attitude_ekf_quaternions_nobias.png

%%
	function x_dot = euler_kinematics(x, pqr)
		e_thta_	= x(2);
		e_phi_	= x(3);
		
		x_dot	= [...
			-sin(e_thta_) 0 1; ...
			sin(e_phi_)*cos(e_thta_) cos(e_phi_) 0; ...
			cos(e_phi_)*cos(e_thta_) -sin(e_phi_) 0] \ pqr;
	end

	function x_dot = quaternion_kinematics(x, pqr)
		p_ = pqr(1);
		q_ = pqr(2);
		r_ = pqr(3);
		x_dot = 0.5*[0 -p_ -q_ -r_; p_ 0 r_ -q_; q_ -r_ 0 p_; r_ q_ -p_ 0]*x;
	end

	function quat_ = euler2quat(euler_angles_)
		e_psi_	= euler_angles_(1);
		e_thta_	= euler_angles_(2);
		e_phi_	= euler_angles_(3);
		quat_	=	[ ...
			cos(e_psi_/2)*cos(e_thta_/2)*cos(e_phi_/2) + sin(e_psi_/2)*sin(e_thta_/2)*sin(e_phi_/2); ...
			cos(e_psi_/2)*cos(e_thta_/2)*sin(e_phi_/2) - sin(e_psi_/2)*sin(e_thta_/2)*cos(e_phi_/2); ...
			cos(e_psi_/2)*sin(e_thta_/2)*cos(e_phi_/2) + sin(e_psi_/2)*cos(e_thta_/2)*sin(e_phi_/2); ...
			sin(e_psi_/2)*cos(e_thta_/2)*cos(e_phi_/2) - cos(e_psi_/2)*sin(e_thta_/2)*sin(e_phi_/2)];
 
	end

	function euler_angles_ = quat2euler(quat_)
		e0_ = quat_(1);		e1_ = quat_(2);
		e2_ = quat_(3);		e3_ = quat_(4);
		e_psi_	= atan2( 2*(e0_*e3_ + e1_*e2_), e0_^2 + e1_^2 - e2_^2 - e3_^2 );
		e_thta_ = asin( 2*(e0_*e2_ - e1_*e3_) );
		e_phi_	= atan2( 2*(e0_*e1_ + e2_*e3_), e0_^2 + e3_^2 - e1_^2 - e2_^2 );
		
		euler_angles_ = [e_psi_; e_thta_; e_phi_];
	end
end