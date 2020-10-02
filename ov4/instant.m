close all;
clear;

%% Insert true system
m = 20;
beta = 0.1;
k = 5;

A = [0 1;
     -k/m -beta/m];              
B = [0; 1/m];              

%% Define filters
lambda_0 = 0.5;
lambda_1 = 0.5;

[ ~, ~,C_f_1,D_f_1] = tf2ss([0 0 1],[1 lambda_1 lambda_0]);
[ ~, ~,C_f_2,D_f_2] = tf2ss([0 1 0],[1 lambda_1 lambda_0]);
[A_f,B_f,C_f_3,D_f_3] = tf2ss([1 0 0],[1 lambda_1 lambda_0]);                 
% A_f and B_f is only needed once as it is the same all over. Hence tildes

%% Adaptive gain
gamma   = eye(3);                              
 
%% Simulation MATLAB
h       = 0.01;  % sample time (s)
N       = 50000; % number of samples

t       = 0:h:h*(N-1);

% Define input as a function of t
u       = sin(t);                               


% Memory allocation
x       = zeros(2, N);
x_z     = zeros(2, 1);
x_phi   = zeros(2, 1);
theta   = zeros(3, N);
U       = zeros(N-1,1);
m_plot = zeros(N, 1); 

% Initial estimates
theta(:,1) = [15 0.5 3]';  


% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1
    
     if t(n) > 20
        m = 20*(2-exp(-0.01*(t(n)-20)));
        m_plot(n) = m;
    else
        m = 20;
        m_plot(n) = m;
    end
    
    A = [0 1;
     -k/m -beta/m];
    B = [0; 1/m];   
 
    
    % Simulate true system
    x_dot           = A*x(:, n) + B*u(n);
    x(:, n+1)       = x(:, n) + h*x_dot;
    y               = x(1, n);
    
    % Generate z and phi by filtering known signals
    x_z_n           = x_z + (A_f*x_z + B_f*u(n))*h;        % u is unfiltered 'z'
    z               = C_f_1*x_z;                           %  1/Lambda * u
    
    x_phi_n         = x_phi + (A_f*x_phi + B_f*y)*h;
                    %      s^2/Lambda * y           s/Lambda * y               1/Lambda * y
    phi             = [(C_f_3*x_phi + D_f_3*y); (C_f_2*x_phi + D_f_2*y); (C_f_1*x_phi + D_f_1*y)];
    
    % Calculate estimation error
    epsilon         = z-theta(:, n)'*phi;  
    
    % Update law
    theta_dot       = gamma*epsilon*phi;                   
    theta(:, n+1)   = theta(:, n) + theta_dot*h;
    
    % Set values for next iteration
    x_phi           = x_phi_n;
    x_z             = x_z_n;
end

% Plots
figure
subplot(3,1,1)
plot(t, theta(1,:)); hold on
plot(t, m_plot); hold off
ylabel('m')
legend('estimate','true value')
grid
subplot(3,1,2)
plot(t, theta(2,:)); hold on
plot([t(1), t(end)],[beta, beta]); hold off
ylabel('\beta')
legend('estimate','true value')
grid
subplot(3,1,3)
plot(t, theta(3,:)); hold on
plot([t(1), t(end)],[k, k]); hold off
ylabel('k')
grid
legend('estimate','true value')
sgtitle('Parameter estimates')
xlabel('t [s]')