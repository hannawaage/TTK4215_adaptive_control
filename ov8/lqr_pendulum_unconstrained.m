close all;
clear;

% True parameters
m    = 0.2;    % [kg]
l    = 1;      % [m]
g    = 9.81;   % [m/sec^2]
beta = 0.003;

% Linearization
A = [0, 1; g/l, -beta/(m*l^2)];
B = [0; 1/(m*l^2)];

Co = ctrb(A,B)
rank(Co)

Q = diag([1,1]);
R = 1;

[K,S,e] = lqr(A,B,Q,R)
pendulum_dyn = @(t,x) [x(2); (1/(m*l^2))*(-K*[x(1)-pi; x(2)] -beta*x(2) - m*g*l*sin(x(1)))];


x_0 = [160*pi/180;30*pi/180];
simTime = 10; % [s]
[t,x] = ode45(pendulum_dyn,[0,simTime],x_0);


control_fun = @(x) -K(1)*(x(:,1)-pi)-K(2)*x(:,2);

figure
subplot(3,1,1)
plot(t,x(:,1)*180/pi)
ylabel('q')
subplot(3,1,2)
plot(t,x(:,2)*180/pi)
ylabel('dq/dt')
subplot(3,1,3)
plot(t,control_fun(x))
ylabel('\tau')
xlabel('Time [s]')




