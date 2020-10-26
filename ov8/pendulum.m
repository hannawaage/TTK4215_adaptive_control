function dx = pendulum(t,x,m,l,beta,g, A,B,C_z,D_z,C_phi_1,C_phi_2,C_phi_3,Gamma)
% True system states
% x(1) = q
% x(2) = dq
% Filter states
% x(3): Internal state for filters 1 and 4 (z and phi_3)
% x(4): Internal state for filters 1 and 4 (z and phi_3)
% x(5): Internal state for filter 2 (phi_1)
% x(6): Internal state for filter 2
% x(7): Internal state for filter 3 (phi_2)
% x(8): Internal state for filter 3
% Estimates
% x(9):  theta_1 (1/(ml^2))
% x(10): theta_2 (g/l)
% x(11): theta_3 (beta/(ml^2))

q = x(1);

% z = s^2/Lambda(s) q 
z = C_z*x(3:4) + D_z*q;

% phi_1 = 1/Lambda(s) tau
phi_1 = C_phi_1*x(5:6);

% phi_2 = -1/Lambda(s) sin(q)
phi_2 = C_phi_2*x(7:8);

% phi_3 = -s/Lambda(s) q
phi_3 = C_phi_3*x(3:4); % NB: share states with z-filter

phi = [phi_1; phi_2; phi_3];
theta = x(9:11);
eps = (z-theta'*phi)/(1+phi'*phi);

l_hat = g/theta(1);
m_hat = theta(2)^2/(theta(1)*g^2);
beta_hat = theta(3)/theta(1);

% Linearization
A_hat = [0, 1; g/l_hat, -beta_hat/(m_hat*l_hat^2)];
B_hat = [0; 1/(m_hat*l_hat^2)];

Q = diag([10,1]);
R = 1;

K_hat = dlqr(A_hat,B_hat,Q,R);

q_dot = x(2);
E_tilde = 0.5*m_hat*l^2*q_dot^2 - m_hat*g*l_hat*(1 + cos(q));
ntimes = floor(abs(q/(2*pi)));
close_enough = (3/4*pi*(1 + ntimes) <= abs(q)) && (abs(q) <= 5/4*pi*(1 + ntimes));
    
if close_enough
    tau = -K_hat*[q-pi; q_dot];
else
    tau = beta*q_dot - E_tilde*q_dot;
end
    
if abs(tau) > 0.5
    tau = sign(tau)*0.5;
end

dx_1 = x(2);
dx_2 = 1/(m*l^2) * tau - (g/l)*sin(q) - beta/(m*l^2)*x(2);

dx_34 = A*x(3:4) + B*q;
dx_56 = A*x(5:6) + B*tau;
dx_78 = A*x(7:8) + B*sin(q);

dtheta = Gamma*eps*phi;

use_projection = true;

% We know all values should be positive.
% Use projection to bound from below
min = 0;


if(use_projection)
    % Assuming diagonal Gamma. Decouple gradient with projection.
    if (theta(1) > min) || (dtheta(1) > 0)
        % In valid region or derivative pointing in towards valid region.
        % No correction needed
    else
        dtheta(1) = 0;
    end
    
    % Assuming diagonal Gamma. Decouple gradient with projection.
    if (theta(2) > min) || (dtheta(2) > 0)
        % In valid region or derivative pointing in towards valid region.
        % No correction needed
    else
        dtheta(2) = 0;
    end
    
    % Assuming diagonal Gamma. Decouple gradient with projection.
    if (theta(3) > min) || (dtheta(3) > 0)
        % In valid region or derivative pointing in towards valid region.
        % No correction needed
    else
        dtheta(3) = 0;
    end
    
    
end

dx = [dx_1; dx_2; dx_34; dx_56; dx_78; dtheta];
end

