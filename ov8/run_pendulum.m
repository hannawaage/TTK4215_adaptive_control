% Assignment 7, Problem 1: Pendulum dynamics
close all;
clear;

% True parameters
m    = 0.2;    % [kg]
l    = 1;   % [m]
g    = 9.81;   % [m/sec^2]
beta = 0.003;

% Initial guesses
m0 = 0.1;
l0 = 1.2;
beta0 = 0.006;

eps = 0.009*pi;
q_0 = pi + eps;
q_dot_0 = 0;

% Initial state
x0 = [q_0 - pi; q_dot_0; 0; 0; 0; 0; 0; 0; 1/(m0*l0^2); g/l0; beta0/(m0*l0^2)];
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

%% Filter design

% Frequency of pendulum:
omega = 1/(sqrt(l/g));

% 1/(s+a)^2
a = 1;
b = 1;

num1 = [1, 0, 0];
num2 = [1];
num3 = [-1];
num4 = [-1 0];
den = [1, a+b, a*b];
[A,B,C_z,D_z] = tf2ss(num1, den)
[~,~,C_phi_1,~] = tf2ss(num2, den);
[~,~,C_phi_2,~] = tf2ss(num3, den);
[~,~,C_phi_3,~] = tf2ss(num4, den);

%figure
%bode(tf(num2,den))

%% Adaptive gains
Gamma = diag([5000, 50, 5]);

%% Simulate
simTime = 50; % [s]
[t,x] = ode45(@(t,y) pendulum(t,y,m,l,beta,g,A,B,C_z,D_z,C_phi_1,C_phi_2,C_phi_3,Gamma),[0,simTime],x0);

%% Plot
figure
subplot(2,1,1)
plot(t,x(:,1) + pi)
ylabel('Angle')
subplot(2,1,2)
plot(t,x(:,1))
ylabel('Angular velocity')

figure
subplot(3,1,1)
plot(t,x(:,9)); hold on
plot([t(1) t(end)],[1/(m*l^2),1/(m*l^2)]); hold off
title('\theta vs \theta^*')
ylabel('\theta_1^* = 1/(ml^2)')
subplot(3,1,2)
plot(t,x(:,10)); hold on
plot([t(1) t(end)],[g/l,g/l]); hold off
ylabel('\theta_{2_1}^* = g/l')
subplot(3,1,3)
plot(t,x(:,11)); hold on
plot([t(1) t(end)],[beta/(m*l^2),beta/(m*l^2)]); hold off
xlabel('Time [s]')
ylabel('\theta_{2_2}^* = \beta/(ml^2)')
