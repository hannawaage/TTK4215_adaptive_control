close all;
clear;

%% Insert true system
% Task i) 20 MPH: 

k_0 = 0.82;
w_0 = 20;
zeta_0 = 0.32;
k_1 = 0.055;
w_1 = 14.2;
zeta_1 = 0.318;


% tf theta_p -> r = G0/(1+G0)
b_0 = k_0*w_0^2;
a_0 = [1 2*zeta_0*w_0 w_0^2];
[A_0, B_0, C_0, D_0] = tf2ss(b_0, a_0);

% tf theta_dot -> theta_p = G1
b_1 = k_1*w_1^2;
a_1 = [1 2*zeta_1*w_1 w_1];
[A_1, B_1, C_1, D_1] = tf2ss(b_1, a_1);

%% Define filters
lambda_0_0 = 100;
lambda_0_1 = 500;

[ ~, ~,C_f_1_0,D_f_1_0] = tf2ss([0 0 1],[1 lambda_0_1 lambda_0_0]);
[ ~, ~,C_f_2_0,D_f_2_0] = tf2ss([0 1 0],[1 lambda_0_1 lambda_0_0]);
[A_f_0,B_f_0,C_f_3_0,D_f_3_0] = tf2ss([1 0 0],[1 lambda_0_1 lambda_0_0]);

lambda_1_0 = 100;
lambda_1_1 = 500;

[ ~, ~,C_f_1_1,D_f_1_1] = tf2ss([0 0 1],[1 lambda_1_1 lambda_1_0]);
[ ~, ~,C_f_2_1,D_f_2_1] = tf2ss([0 1 0],[1 lambda_1_1 lambda_1_0]);
[A_f_1,B_f_1,C_f_3_1,D_f_3_1] = tf2ss([1 0 0],[1 lambda_1_1 lambda_1_0]);


%% Adaptive gain
v_0 = [75 50 100];
Gamma_0   = diag(v_0);        
v_1 = [250 200 100];
Gamma_1 = diag(v_1);

k_max = 1;
w_max = 30;

%% Simulate system
h = 0.001; % sample time (s)
N = 50000; % number of samples
t = 0:h:h*(N-1);

% Define input as a function of t
r = 10*sin(0.2*t) + 0.1396;

% Memory allocation
x = zeros(2, N); 
gamma_0 = zeros(3, N);
gamma_1 = zeros(3, N);
estimates_0 = zeros(3, N);
estimates_1 = zeros(3, N);
x_0 = zeros(2, N);
x_1 = zeros(2, N);
x_z_0 = zeros(2, 1);
x_phi_0 = zeros(2, 1);
x_z_1 = zeros(2, 1);
x_phi_1 = zeros(2, 1);


% Initial estimate
k_0_i = 0.3; w_0_i = 14; zeta_0_i = 0.1;
gamma_0(:,1) = [1/(k_0_i*w_0_i^2); 2*zeta_0_i/(k_0_i*w_0_i); 1/k_0_i];
k_1_i = 0.1; w_1_i = 10; zeta_1_i = 0.1;
gamma_1(:,1) = [1/(k_1_i*w_1_i^2); 2*zeta_1_i/(k_1_i*w_1_i); 1/k_1_i];

% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1

    % Simulate true system
    
    % tf theta_p -> r = G0/(1+G0)
    x_0_dot = A_0*x_0(:, n) + B_0*r(n);
    x_0(:, n+1) = x_0(:, n) + h*x_0_dot;
    y_0 = C_0*x_0(:, n+1) + D_0*r(n);
    theta_p = y_0;
    
    % tf theta_dot -> theta_p = G1
    x_1_dot = A_1*x_1(:, n) + B_1*theta_p;
    x_1(:, n+1) = x_1(:, n) + h*x_1_dot;
    y_1 = C_1*x_1(:, n+1) + D_1*theta_p;
    theta_dot = y_1;
    
    % Generate z_0, z_1, phi_0 and phi_1 by filtering known signals
    x_z_0_n = x_z_0 + (A_f_0*x_z_0 + B_f_0*r(n))*h;
    z_0 = C_f_1_0*x_z_0;
    
    x_phi_0_n = x_phi_0 + (A_f_0*x_phi_0 + B_f_0*theta_p)*h;
    phi_0 = [(C_f_3_0*x_phi_0 + D_f_3_0*theta_p);           %s^2/Lambda *theta_p
            (C_f_2_0*x_phi_0 + D_f_2_0*theta_p);            %s/Lambda *theta_p
            (C_f_1_0*x_phi_0 + D_f_1_0*theta_p)];           %1/Lambda *theta_p
    
    x_z_1_n = x_z_1 + (A_f_1*x_z_1 + B_f_1*theta_p)*h;
    z_1 = C_f_1_1*x_z_1;
    
    x_phi_1_n = x_phi_1 + (A_f_1*x_phi_1 + B_f_1*theta_dot)*h;
    phi_1 = [(C_f_3_1*x_phi_1 + D_f_3_1*theta_dot);
            (C_f_2_1*x_phi_1 + D_f_2_1*theta_dot);
            (C_f_1_1*x_phi_1 + D_f_1_1*theta_dot)];
    
    % Calculate estimation error
    n_s_0 = phi_0'*phi_0;
    n_s_1 = phi_1'*phi_1;
    epsilon_0 = (z_0-gamma_0(:, n)'*phi_0)/(1+n_s_0^2);
    epsilon_1 = (z_1-gamma_1(:, n)'*phi_1)/(1+n_s_1^2);
    
    %Find value of g and grad_g
    k_0_n = 1/gamma_0(3, n);
    w_0_n = sqrt(gamma_0(3, n)/gamma_0(1, n));
    zeta_0_n = 0.5*(gamma_0(2, n)*k_0_n*w_0_n);
    estimates_0(:, n) = [k_0_n, w_0_n, zeta_0_n]';
    
    k_1_n = 1/gamma_1(3, n);
    w_1_n = sqrt(gamma_1(3, n)/gamma_1(1, n));
    zeta_1_n = 0.5*(gamma_1(2, n)*k_1_n*w_1_n);
    estimates_1(:, n) = [k_1_n, w_1_n, zeta_1_n]';
    
    args_0 = [1/k_max - gamma_0(3, n); 1/(k_max*w_max^2) - gamma_0(1, n)];
    args_1 = [1/k_max - gamma_1(3, n); 1/(k_max*w_max^2) - gamma_1(1, n)];
    
    g_0 = max(args_0);
    g_1 = max(args_1);
    grad_g_0 = [-1 0 -1]';
    grad_g_1 = [-1 0 -1]';
    
    if (1/k_max == gamma_0(3, n))
        grad_g_0(1) = 0;
    end
    if (1/(k_max*w_max^2) == gamma_0(1, n))
        grad_g_0(3) = 0;
    end
    
    if (1/k_max == gamma_1(3, n))
        grad_g_1(1) = 0;
    end
    if (1/(k_max*w_max^2) == gamma_1(1, n))
        grad_g_1(3) = 0;
    end
    
    % Update law: gamma0_dot
    if g_0 < 0 || ((g_0 == 0) && ((Gamma_0*epsilon_0*phi_0)'*grad_g_0 <= 0))
        gamma_0_dot = Gamma_0*epsilon_0*phi_0;
        gamma_0(:, n+1) = gamma_0(:, n) + gamma_0_dot*h;
    else
        gamma_0_dot = Gamma_0*epsilon_0*phi_0 - Gamma_0*(grad_g_0'*Gamma_0*grad_g_0)\(grad_g_0*grad_g_0')*Gamma_0*epsilon_0*phi_0;
        gamma_0(:, n+1) = gamma_0(:, n) + gamma_0_dot*h;
    end
    
    % Update law: gamma1_dot
    if g_1 < 0 || ((g_1 == 0) && ((Gamma_1*epsilon_1*phi_1)'*grad_g_1 <= 0))
        gamma_1_dot = Gamma_1*epsilon_1*phi_1;
        gamma_1(:, n+1) = gamma_1(:, n) + gamma_1_dot*h;
    else
        gamma_1_dot = Gamma_1*epsilon_1*phi_1 - Gamma_1*(grad_g_1'*Gamma_1*grad_g_1)\(grad_g_1*grad_g_1')*Gamma_1*epsilon_1*phi_1;
        gamma_1(:, n+1) = gamma_1(:, n) + gamma_1_dot*h;
    end
    
    % Set values for next iteration
    x_z_0 = x_z_0_n;
    x_phi_0 = x_phi_0_n;
    
    x_z_1 = x_z_1_n;
    x_phi_1 = x_phi_1_n;
end


% Plots
figure(1);
sgtitle('Estimates of k_0, \omega_0, \zeta_0');

subplot(3,1,1)
plot(t, estimates_0(1,:)); hold on
plot([t(1), t(end)],[k_0, k_0]); hold off
ylabel('k_0')
grid
legend('estimate','true value')
xlabel('t [s]')

subplot(3,1,2)
plot(t, estimates_0(2,:)); hold on
plot([t(1), t(end)], [w_0 w_0]); hold off
ylabel('\omega_0')
legend('estimate','true value')
grid

subplot(3,1,3)
plot(t, estimates_0(3,:)); hold on
plot([t(1), t(end)],[zeta_0, zeta_0]); hold off
ylabel('\zeta_0')
legend('estimate','true value')
grid

figure(2);
sgtitle('Estimates of k_1, \omega_1, \zeta_1');

subplot(3,1,1)
plot(t, estimates_1(1,:)); hold on
plot([t(1), t(end)],[k_1, k_1]); hold off
ylabel('k_1')
grid
legend('estimate','true value')
xlabel('t [s]')

subplot(3,1,2)
plot(t, estimates_1(2,:)); hold on
plot([t(1), t(end)], [w_1 w_1]); hold off
ylabel('\omega_1')
legend('estimate','true value')
grid

subplot(3,1,3)
plot(t, estimates_1(3,:)); hold on
plot([t(1), t(end)],[zeta_1, zeta_1]); hold off
ylabel('\zeta_1')
legend('estimate','true value')
grid

