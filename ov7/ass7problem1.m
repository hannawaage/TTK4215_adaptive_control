%close all;
clear;

% True parameters
m       = 0.2;  % [kg]
l       = 1;    % [m]
beta    = 0.003;% [kg m^2/sec]
g       = 9.81; % [m/sec^2]

% Initial state
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1; 4; 1];

% True system states
% x(1) = q
% x(2) = q_dot
% x(3) to x(10) = filter states
% Estimates
% x(11): theta_1 (beta/ml^2)
% x(12): theta_21 (g/l)
% x(13): theta_22 (1/ml^2)


%% Define filters
% l/(s+a)^2
a=7;
num1 = [0, -1, 0];
num2 = [0, 0, -1];
num3 = [0, 0, 1];
num4 = [1, 0, 0];
den  = [1, 2*a, a^2];

[~,~,C_phi_1,~] = tf2ss(num1, den);
[~,~,C_phi_2,~] = tf2ss(num2, den);
[~,~,C_phi_3,~] = tf2ss(num3, den);
[A,B,C_z,D_z]   = tf2ss(num4, den);                 
 
%% Simulation MATLAB
h       = 0.01;   % sample time (s)
N       = 50000; % number of samples
t       = 0:h:h*(N-1);

% Define Gamma etc
Gamma_b = diag([15 15000 25]);
Gamma_d = diag([15 15000 25]);  
A_sys = [0 1; m*g/l -beta/(m*l^2)];
B_sys = [0; 1/(m*l^2)];
Q = diag([5 1]);
R = 1;
K = lqr(A_sys, B_sys, Q,R);

simTime = 50; % [s]
[t,x] = ode45(@(t,y) problem1_b(t,y,m,beta,g,l,A,B,C_phi_1,C_phi_2,C_phi_3,C_z,D_z,Gamma_b),[0,simTime],x0);
%[t,x] = ode45(@(t,y) problem1_d(t,y,m,beta,g,l,A,B,C_phi_1,C_phi_2,C_phi_3,C_z,D_z,Gamma_d, K),[0,simTime],x0);

% Plots
figure('Name', 'Parameter estimates')
subplot(3,1,1)
plot(t,x(:,11)); hold on
plot([t(1) t(end)],[beta/m*l^2,beta/m*l^2]); hold off
title('Result of parameter estimation')
ylabel('\theta_1^* = \beta/m*l^2')
subplot(3,1,2)
plot(t,x(:,12)); hold on
plot([t(1) t(end)],[g/l,g/l]); hold off
ylabel('\theta_2^* = g/l')
subplot(3,1,3)
plot(t,x(:,13)); hold on
plot([t(1) t(end)],[1/m*l^2,1/m*l^2]); hold off
xlabel('Time [s]')
ylabel('\theta_3^* = 1/m*l^2')

figure(2)
subplot(2,1,1)
plot(t,x(:,1)); hold on
plot([t(1) t(end)],[pi,pi]); hold off
ylabel('q') 
xlabel('Time [s]')
subplot(2,1,2)
plot(t,x(:,2)); hold on
plot([t(1) t(end)],[0,0]); hold off
ylabel('q_{dot}') 
xlabel('Time [s]')
