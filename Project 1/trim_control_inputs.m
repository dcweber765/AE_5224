function tci = trim_control_inputs(alpha_star,beta_star,g,V_a_star,rho_star,J_x,J_z,J_xz,c,m,b,S,k_motor,S_prop,C_prop,v_star,w_star,C_m_0, C_m_alpha,C_m_q,C_m_delta_e,C_p_r,C_r_r,C_X,C_X_q,C_X_del_e,theta_star,p_star,q_star,r_star,C_p_delta_a,C_p_delta_r, C_r_delta_a,C_r_delta_r,C_p_0,C_r_0,C_p_p,C_r_p,C_p_beta,C_r_beta,Gamma_1,Gamma_2,Gamma_7)
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
tci(1) = delta_e_star;
tci(2) = delta_p_star;
tci(3) = delta_a_star;
tci(4) = delta_r_star;
end