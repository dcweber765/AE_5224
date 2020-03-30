close all; clear variables

x0 = [0;0;100;0;0;.1;0;0;0;0;0;0];

n_pts	= 1001;
time_pts= linspace(0, 300, 1001);
Xstate_from_true	= zeros(12, n_pts);
control_in_OMEGA = zeros(4,n_pts);
X_true =	x0;
Xstate_from_true(:, 1)= X_true;
starT1 = 0;%[V_star,alpha_star,beta_star,delta_e_star,delta_t_star,delta_a_star,delta_r_star,u_star,v_star,w_star,p_star,q_star,r_star,psi_star,theta_star,phi_star,m,J_x,J_y,J_z,J_zx,J_xz, rho_star,S,c_bar,b];
%% RK4
for m1 = 2:n_pts
    
    % RK4 step (Euler)
    dt	= time_pts(m1) - time_pts(m1 - 1);
    u1	= control_in_OMEGA(:, m1 - 1);
    u2	= control_in_OMEGA(:, m1);
    u12	= 0.5*(u1 + u2);							% Approx (p,q,r) at (t + 0.5dt)
    
    k1	= dt*quadDynamics(X_true,				u1, starT1);
    k2	= dt*quadDynamics((X_true + 0.5*k1),	u12, starT1);
    k3	= dt*quadDynamics((X_true + 0.5*k2),	u12, starT1);
    k4	= dt*quadDynamics((X_true + k3),		u2, starT1);
    X_true	= X_true + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4;
    
    % Record
    Xstate_from_true(:, m1)	= X_true;
end

figure('Name', 'Pos Stick-Fixed Response RK4')
subplot(311)
plot(time_pts, Xstate_from_true(1,:))
xlabel('t (s)'); ylabel('\Delta p_x m'); grid on;

subplot(312)
plot(time_pts, Xstate_from_true(2,:))
xlabel('t (s)'); ylabel('\Delta p_y m'); grid on;

subplot(313)
plot(time_pts, Xstate_from_true(3,:))
xlabel('t (s)'); ylabel('\Delta p_z m'); grid on;

figure('Name', 'Vel Stick-Fixed Response RK4')
subplot(311)
plot(time_pts, Xstate_from_true(4,:))
xlabel('t (s)'); ylabel('\Delta u m/s'); grid on;

subplot(312)
plot(time_pts, Xstate_from_true(5,:))
xlabel('t (s)'); ylabel('\Delta v m/s'); grid on;

subplot(313)
plot(time_pts, Xstate_from_true(6,:))
xlabel('t (s)'); ylabel('\Delta w m/s'); grid on;


