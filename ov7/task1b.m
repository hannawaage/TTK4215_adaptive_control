close all;
clear;

%% Insert true system
l = 1;
m = 0.2;
beta = 0.003;
grav = 9.81;

theta_real = [1/(m*l^2); beta/(m*l^2); grav/l];

A = [0 1;
     0 -beta/(m*l^2)];              
B = [0 0;
     1/(m*l^2) -grav/l];            

%% Define valid area
l_max = 2;
m_max = 1;

theta_1_min = 1/(m_max*l_max^2);
theta_3_min = grav/(l_max);

%% Define filters
lambda_0 = 25;
lambda_1 = 10;

[ ~, ~,C_f_1,D_f_1] =           tf2ss([0 0 1],[1 lambda_1 lambda_0]);
[ ~, ~, C_f_1_neg, D_f_1_neg] = tf2ss([0 0 -1],[1 lambda_1 lambda_0]);
[ ~, ~,C_f_2,D_f_2] =           tf2ss([0 -1 0],[1 lambda_1 lambda_0]);
[A_f,B_f,C_f_3,D_f_3] =         tf2ss([1 0 0],[1 lambda_1 lambda_0]);                 
% A_f and B_f is only needed once as it is the same all over. Hence tildes

%% Adaptive gain & forgetting factor 
p0 = 1000;
P_0 = p0*eye(3);

zeta = 1.2;

%% Simulate system

h = 0.001; % sample time (s)
N = 50000; % number of samples
t = 0:h:h*(N-1);

% Define input as a function of t
tau = sin(0.2*t) + sin(10*t);

% Memory allocation
x = zeros(2, N); 
theta = zeros(3, N);
estimates = zeros(3, N);
x_z = zeros(2, 1);
x_phi = zeros(3, 1);
x_phi_tau = zeros(2, 1);
x_phi_q = zeros(2, 1);
x_phi_sinq = zeros(2, 1);
P = zeros(3, 3, N);


% Initial estimate
l_0 = 0.3;
m_0 = 0.8;
beta_0 = 0.5;
theta(:, 1) = [1/(m_0*l_0^2); beta_0/(m_0*l_0^2); grav/l_0];
P(:, :, 1) = P_0;

% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N
    
    % Simulate true system
    q = x(1, n);
    u = [tau(n); sin(q)];
    x_dot = A*x(:, n) + B*u;
    x(:, n+1) = x(:, n) + h*x_dot;
    
    % Generate z and phi by filtering known signals
    
    x_z_n = x_z + (A_f*x_z + B_f*q)*h;
    z = C_f_3*x_z_n;
    
    x_phi_tau_n = x_phi_tau + (A_f*x_phi_tau + B_f*tau(n))*h;
    x_phi_q_n = x_phi_q + (A_f*x_phi_q + B_f*q)*h;
    x_phi_sinq_n = x_phi_sinq + (A_f*x_phi_sinq + B_f*sin(q))*h;
    
    phi = [(C_f_1*x_phi_tau_n + D_f_1*tau(n));      %1/Lambda *tau
           (C_f_2*x_phi_q_n + D_f_2*q);              %-s/Lambda *q
           (C_f_1_neg*x_phi_sinq_n + D_f_1_neg*sin(q))];      %-1/Lambda *sin(q)


    % Calculate estimation error
    n_s = phi'*phi;
    m_squared = 1 + n_s^2;
    epsilon = (z-theta(:, n)'*phi)/m_squared;
    
    %Find value of g and grad_g
    estimates(:, n) = [grav/theta(3, n), theta(3, n)^2/(theta(1, n)*grav^2), theta(2, n)/theta(1, n)]; % l, m, beta
    args_test = [estimates(1, n) - l_max, estimates(2, n)- m_max];
    args = [theta(1, n) - theta_1_min , theta(3, n) - theta_3_min ];
    g = max(args_test);
    grad_g = [-1 0 -1]';
    if (args_test(1) == 0)
        grad_g(1) = 0;
    end
    if (args_test(2) == 0)
        grad_g(3) = 0;
    end
    
    % Update law
    if g < 0 || ((g == 0) && (P(:, :, n)*epsilon*phi)'*grad_g <= 0)
        theta_dot = P(:, :, n)*epsilon*phi;
        theta(:, n+1) = theta(:, n) + theta_dot*h;
        P_dot = zeta*P(:, :, n) - P(:, :, n)*(phi*phi')/m_squared*P(:, :, n);
        P(:, :, n+1) = P(:, :, n) + h*P_dot;
    else
        theta_dot = P(:, :, n)*epsilon*phi - P(:, :, n)*(grad_g'*P(:, :, n)*grad_g)\(grad_g*grad_g')*P(:, :, n)*epsilon*phi;
        theta(:, n+1) = theta(:, n) + theta_dot*h;
        P_dot = 0;
        P(:, :, n+1) = P(:, :, n);
    end
    
    % Set values for next iteration
    x_phi_tau = x_phi_tau_n;
    x_phi_q = x_phi_q_n;
    x_phi_sinq = x_phi_sinq_n;
 
    x_z = x_z_n;
end


% Plots
figure

subplot(3,1,1)
plot(t, estimates(1, :)); hold on
plot([t(1), t(end)],[l, l]); hold off
ylabel('l')
grid
legend('estimate','true value')
sgtitle('Parameter estimates')
xlabel('t [s]')

subplot(3,1,2)
plot(t, estimates(2,:)); hold on
plot([t(1), t(end)], [m m]); hold off
ylabel('m')
legend('estimate','true value')
grid

subplot(3,1,3)
plot(t, estimates(3,:)); hold on
plot([t(1), t(end)],[beta, beta]); hold off
ylabel('\beta')
legend('estimate','true value')
grid

