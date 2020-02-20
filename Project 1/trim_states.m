function ts = trim_states(V_a_star,R_star,gamma_star,alpha_star,beta_star,phi_star) 
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

ts(1) = u_star;
ts(2) = v_star;
ts(3) = w_star;
ts(4) = theta_star;
ts(5) = p_star;
ts(6) = q_star;
ts(7) = r_star;
end
