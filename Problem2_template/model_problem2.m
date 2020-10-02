function dxdt = model_problem2(t,x,A,B,gamma_1,gamma_2)
% 1-2: states
% 3-4: state predictions
% 5-8: (flattened, a11, a21, a12, a22) estimates of A
% 9-10: estimates of B

dxdt = zeros(10,1);


u = 10*sin(2*t) + 7*cos(3.6*t);                                %%% Define your input here

% Simulate true system
dxdt(1:2,1) = A*x(1:2,1) + B*u;

% Reshape estimates of A from vector to matrix
A_hat = reshape(x(5:8,1),[2 2]);
B_hat = x(9:10,1);

% Update law for state prediction
x_hat = x(3:4);
dxdt(3:4,1) = A_hat*x_hat + B_hat*u;                      %%% Insert update law for x_hat

% Update law for estimates of A and B
eps_1 = x(1:2,1) - x_hat;
dA_hat = gamma_1*eps_1*x_hat';                           %%% Insert update law for A_hat
dxdt(5:8,1) = dA_hat(:);

dxdt(9:10,1) = gamma_2*eps_1*u';                     %%% Insert update law for B_hat

end

