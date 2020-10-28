close all;
clear;

%% Insert true system
M = 10;
f = 1;
k = 9;

%% Find true system matrices for simulation

[A_p,B_p,C_p,D_p] = tf2ss([0 0 1/M],[1 f/M k/M]);        
[A_m,B_m,C_m,D_m] = tf2ss([0 0 1],[1 sqrt(2) 1]);       

%% Adaptive gains

p0 = 1;

gamma = 100*[1, 1, 1, 1];
Gamma = diag(gamma);

%% Simulate system

h = 0.001; % sample time (s)
N = 50000; % number of samples
t = 0:h:h*(N-1);

% Define reference
%r = 10;
r = 2*sin(3*t) + 5*sin(t);

% Memory allocation
x_p = zeros(2, N);
x_m = zeros(2, N);
theta = zeros(4, N);
w = zeros(4, N);
phi = zeros(4, N);
y_p_t = zeros(1, N);
y_m_t = zeros(1, N);
e_1_t = zeros(1, N);


% Initial estimate
x_p(:, 1) = 0;
x_m(:, 1) = 0;
theta(:, 1) = 0;
w(:, 1) = 0;
phi(:, 1) = 0;
u_p = 0;

% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1
    
    x_p_dot = A_p*x_p(:, n) + B_p*u_p;
    x_p(:, n+1) = x_p(:, n) + h*x_p_dot;
    y_p = C_p*x_p(:, n+1);
    y_p_t(n) = y_p;
    
    x_m_dot = A_m*x_m(:, n) + B_m*r(n);
    x_m(:, n+1) = x_m(:, n) + h*x_m_dot;
    y_m = C_m*x_m(:, n+1);
    y_m_t(n) = y_m;
    
    e_1 = y_p - y_m;
    e_1_t(n) = e_1;
    
    w_1_dot = -w(1, n) + u_p;
    w_2_dot = -w(2, n) + y_p;
    
    w_1 = w(1, n) + h*w_1_dot;
    w_2 = w(2, n) + h*w_2_dot;
    w(:, n + 1) = [w_1, w_2, y_p, r(n)];
    
    phi_dot = -p0*phi(:, n) + w(:, n);
    phi(:, n+1) = phi(:, n) + h*phi_dot;
    
    u_p = theta(:, n)'*w(:, n) - phi(:, n)'*Gamma*phi(:, n)*e_1;
    theta_dot = -Gamma*e_1*phi(:, n);
    theta(:, n + 1) = theta(:, n) + h*theta_dot;
    
    
end


figure(2);
tit = 'System trajectories with r = 2*sin(3*t) + 5*sin(t)';
sgtitle(tit);
subplot(2, 1, 1);
plot(t, y_p_t); hold on;
plot(t, y_m_t); hold on;
plot(t, r); hold off
ylabel('y')
legend('y_p(t)','y_m(t)', 'r')
grid

subplot(2, 1, 2);
plot(t, e_1_t); hold on;
ylabel('e_1');
legend('e_1');
grid


