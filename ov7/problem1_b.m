function dx = problem1_b(t, x, m, beta, g, l, A, B, C_phi_1, C_phi_2, C_phi_3, C_z, D_z, Gamma)
% True system states
% x(1) = q
% x(2) = q_dot
% Filter states
% x(3): Internal state for filter 1
% x(4): Internal state for filter 1
% x(5): Internal state for filter 2
% x(6): Internal state for filter 2
% x(7): Internal state for filter 3
% x(8): Internal state for filter 3
% x(9): Internal state for filter 4
% x(10): Internal state for filter 4
% Estimates
% x(11): theta_1 (beta/ml^2)
% x(12): theta_21 (g/l)
% x(13): theta_22 (1/ml^2)
q = x(1);
tau = 8*sin(1.34*t)+10*sin(3.12*t) + sin(2*t);
% Parameter estimator for b/ml^2, g/l and 1/ml^2
% theta_1^* = [beta/ml^2, g/l, 1/ml^2]^T
% z = s^2/Lambda(s) q
z = C_z*x(9:10) + D_z*x(1);
% phi_1 = -s/Lambda(s) q
phi_1 = C_phi_1*x(3:4);
% phi_2 = -1/Lambda(s) sin(q)
phi_2 = C_phi_2*x(5:6);
% phi_3 = 1/Lambda(s) tau
phi_3 = C_phi_3*x(7:8);
phi = [phi_1; phi_2; phi_3];
theta = x(11:13);
epsilon = (z-theta'*phi)/(1+phi'*phi);
dx_1 = x(2);
dx_2 = tau/(m*l^2)-beta/(m*l^2)*x(2)-(g/l)*sin(x(1));
dx_34 = A*x(3:4) + B*q;
dx_56 = A*x(5:6) + B*sin(q);
dx_78 = A*x(7:8) + B*tau;
dx_910 = A*x(9:10) + B*q;
dtheta = Gamma*epsilon*phi;
dx = [dx_1; dx_2; dx_34; dx_56; dx_78; dx_910; dtheta];
end
