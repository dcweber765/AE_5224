function [A_lat,B_lat] = lateral_linear_model(Yv,Yp,Yr,Lv,Lp,Lr,Nv,Np,Nr,Y_delta_a,Y_delta_r,L_delta_a,L_delta_r,N_delta_a,N_delta_r,theta_star,phi_star,p_star,q_star,r_star,g)
Lat1A = Yv;
Lat2A = Yp;
Lat3A = Yr;
Lat4A = g*cos(theta_star)*cos(phi_star);
Lat1_5A = 0;
Lat5A = Lv;
Lat6A = Lp;
Lat7A = Lr;
Lat8A = 0;
Lat2_5A = 0;
Lat9A = Nv;
Lat10A = Np;
Lat11A = Nr;
Lat12A = 0;
Lat3_5A = 0;
Lat13A = 0;
Lat14A = 1;
Lat15A = cos(phi_star)*tan(theta_star);
Lat16A = q_star*cos(phi_star)*tan(theta_star) - r_star*sin(phi_star)*tan(theta_star);
Lat4_5A = 0;
Lat5_1A = 0;
Lat5_2A = 0;
Lat5_3A = cos(phi_star)*1/cos(theta_star);
Lat5_4A = p_star*cos(phi_star)*1/cos(theta_star) - r_star*sin(phi_star)*1/cos(theta_star);
Lat5_5A = 0;

A_lat = [Lat1A, Lat2A, Lat3A, Lat4A Lat1_5A;
          Lat5A, Lat6A, Lat7A, Lat8A Lat2_5A;
          Lat9A, Lat10A, Lat11A, Lat12A Lat3_5A;
          Lat13A, Lat14A, Lat15A, Lat16A Lat4_5A;
          Lat5_1A Lat5_2A Lat5_3A Lat5_4A Lat5_5A]
      
B_lat = [Y_delta_a Y_delta_r;
         L_delta_a  L_delta_r;
         N_delta_a  N_delta_r;
         0 0;
         0 0]
end