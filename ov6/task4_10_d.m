close all;
clear;


%% Insert true system
m = 15;
beta = 0.2;
k = 2;

A = [0 0 0;
     0 0 1;
     0 0 -beta/m];              
B = [0; 0; 1/m];              

%% Define filters
lambda_0 = 6;
lambda_1 = 1;

[ ~, ~,C_f_1,D_f_1] = tf2ss([0 0 1],[1 lambda_1 lambda_0]);
[ ~, ~,C_f_2,D_f_2] = tf2ss([0 1 0],[1 lambda_1 lambda_0]);
[A_f,B_f,C_f_3,D_f_3] = tf2ss([1 0 0],[1 lambda_1 lambda_0]);                 
% A_f and B_f is only needed once as it is the same all over. Hence tildes

%% Adaptive gain
gamma_1   = 300;        
v = [500 5];
gamma_2 = diag(v);
%% Simulate system

h = 0.01; % sample time (s)
N = 50000; % number of samples
t = 0:h:h*(N-1);

% Define input as a function of t
u = 5*sin(2*t) + 10.5;

% Memory allocation
x = zeros(3, N); 
theta_1 = zeros(1, N);
x_z_2 = zeros(2, 1);
x_phi_2 = zeros(2, 1);
theta_2 = zeros(2, N);

% Initial estimate
theta_1(:,1) = 1; %k
theta_2(:,1) = [10 ; 1]; %m beta

% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1
    
    % Simulate true system
    x_dot = A*x(:, n) + B*u(n);
    x(:, n+1) = x(:, n) + h*x_dot;
    y_1 = x(1, n);
    y_2 = x(2, n);
    dy_2 = x(3, n);
    
    % Generate z_2 and phi_2 by filtering known signals
    x_z_2_n = x_z_2 + (A_f*x_z_2 + B_f*u(n))*h;
    z_2 = C_f_1*x_z_2;
    x_phi_2_n = x_phi_2 + (A_f*x_phi_2 + B_f*y_2)*h;
    phi_2 = [(C_f_3*x_phi_2 + D_f_3*y_2);
            (C_f_2*x_phi_2 + D_f_2*y_2)];

    % Generate z_1 and phi_1
    z_1 = u(n);
    phi_1 = y_1 - y_2;
    
    % Calculate estimation error
    n_s_1 = phi_1'*phi_1;
    n_s_2 = phi_2'*phi_2;
    epsilon_1 = (z_1-theta_1(:, n)'*phi_1)/(1+n_s_1^2);
    epsilon_2 = (z_2-theta_2(:, n)'*phi_2)/(1+n_s_2^2);
    
    %Find value of g and grad_g
    k_n = theta_1(:, n);
    m_n = theta_2(1, n);
    beta_n = theta_2(2, n);
    args = [-beta_n, beta_n - 1, 0.1 - k_n, 10 - m_n];
    g = max(args);
    grad_g_2 = [-1 -1]';
    if (k_n == 0.1)
        grad_g_2(1) = 0;
    end
    if (m_n == 10)
        grad_g_2(2) = 0;
    end
    
    % Update law
    if g < 0 || ((g == 0) && ((gamma_2*epsilon_2*phi_2)'*grad_g_2 <= 0))
        theta_2_dot = gamma_2*epsilon_2*phi_2;
        theta_2(:, n+1) = theta_2(:, n) + theta_2_dot*h;
        theta_1_dot = gamma_1*epsilon_1*phi_1;
        theta_1(:, n+1) = theta_1(:, n) + theta_1_dot*h;
    else
        theta_1_dot = gamma_1*epsilon_1*phi_1;
        theta_1(:, n+1) = theta_1(:, n) + theta_1_dot*h;
        
        theta_2_dot = gamma_2*epsilon_2*phi_2 - gamma_2*(grad_g_2'*gamma_2*grad_g_2)\(grad_g_2*grad_g_2')*gamma_2*epsilon_2*phi_2;
        theta_2(:, n+1) = theta_2(:, n) + theta_2_dot*h;
    end
    
    
    % Set values for next iteration
    x_phi_2 = x_phi_2_n;
    x_z_2 = x_z_2_n;
end


% Plots
figure

subplot(3,1,3)
plot(t, theta_1(1,:)); hold on
plot([t(1), t(end)],[k, k]); hold off
ylabel('k')
grid
legend('estimate','true value')
sgtitle('Parameter estimates')
xlabel('t [s]')

subplot(3,1,1)
plot(t, theta_2(1,:)); hold on
plot([t(1), t(end)], [m m]); hold off
ylabel('m')
legend('estimate','true value')
grid

subplot(3,1,2)
plot(t, theta_2(2,:)); hold on
plot([t(1), t(end)],[beta, beta]); hold off
ylabel('\beta')
legend('estimate','true value')
grid

