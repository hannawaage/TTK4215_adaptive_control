close all;
clear;

%% Insert true system
l = 1;
m = 0.2;
beta = 0.003;
g = 9.81;

%% Find LQR gain

A_lin = [0 1;
        g/l -beta/(m*l^2)];              
B_lin = [0 ;
        1/(m*l^2)];
    
Q = eye(2);
Q(1) = 10;
R = 0.1;
K = lqr(A_lin, B_lin, Q, R);

%% Simulate system

h = 0.001; % sample time (s)
N = 50000; % number of samples
t = 0:h:h*(N-1);

% Memory allocation
x = zeros(2, N); 
u_t = zeros(1, N);

% Initial estimate
eps = 0.7*pi;
q_0 = pi + eps;
q_dot_0 = -5*eps;
x_0 = [q_0 - pi, q_dot_0];
x(:, 1) = x_0;

% Main loop. Simulate using forward Euler (x[k+1] = x[k] + h*x_dot)
for n = 1:N-1
    q = x(1, n) + pi;
    q_dot = x(2, n);
    
    %% a)
    %u = -K*[q-pi; q_dot];
    
    %% b)
%     u = -K*[q-pi; q_dot];
%     if abs(u) > 0.5
%         u = sign(u)*0.5;
%     end

    %% e)
%     E_tilde = 0.5*m*l^2*q_dot^2 - m*g*l*(1 + cos(q));
%     u = beta*q_dot - E_tilde*q_dot;
%     if abs(u) > 0.5
%         u = sign(u)*0.5;
%     end
    
    %% f)
    E_tilde = 0.5*m*l^2*q_dot^2 - m*g*l*(1 + cos(q));
    ntimes = floor(abs(q/(2*pi)));
    close_enough = (3/4*pi*(1 + ntimes) <= abs(q)) && (abs(q) <= 5/4*pi*(1 + ntimes));
    
    if close_enough
        u = -K*[q-pi; q_dot];
    else
        u = beta*q_dot - E_tilde*q_dot;
    end
    
    if abs(u) > 0.5
        u = sign(u)*0.5;
    end
    
    u_t(n) = u;
    
    dx_1 = q_dot;
    dx_2 = u/(m*l^2)-beta/(m*l^2)*q_dot-(g/l)*sin(q);
    x_dot = [dx_1; dx_2];
    x(:, n+1) = x(:, n) + h*x_dot;
end


figure(2);
tit = 'System trajectories with epsilon = ' + string(eps/pi) + '\pi => q0 = ' + string(q_0/pi) + '\pi and qdot_0 = ' + string(x_0(2));
sgtitle(tit);
subplot(3, 1, 1);
plot(t, x(1, :) + pi); hold on
plot([t(1), t(end)],[pi, pi]); hold off
ylabel('q(t)')
legend('q(t)','reference')
grid

subplot(3, 1, 2);
plot(t, x(2, :)); hold on
plot([t(1), t(end)],[0, 0]); hold off
ylabel('qdot(t)')
legend('qdot(t)','reference');
grid

subplot(3, 1, 3);
plot(t, u_t); hold on
ylabel('u(t)')
legend('u(t)');
grid

