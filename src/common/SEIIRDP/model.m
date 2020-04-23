function [S, E, Ia, Iq, R, D, P] = model(alpha, beta, gamma, delta, lambda0, kappa0, Npop, E0, Ia0, Iq0, R0, D0, t)
% Simulate the time-histories of an epidemic outbreak using generalized SEIIRD model

%% Initial conditions
N = numel(t);  % # time steps
Y = zeros(7, N);  % 6 states (S, E, Ia, Iq, R, D)
Y(1,1) = Npop - E0 - Ia0 - Iq0 - R0 - D0;  % susceptible
Y(2,1) = E0;  % exposed
Y(3,1) = Ia0;  % infectious asymptomatic
Y(4,1) = Iq0;  % quarantined (confirmed) = quarantined + hospitalized
Y(5,1) = R0;  % recovered
Y(6,1) = D0;  % dead
Y(7,1) = 0;  % not susceptible

%% Matrix version of gSEIR model
[Y] = simulate(alpha, beta, gamma, delta, lambda0, kappa0, Y, Npop, t);

S = Y(1, 1:N);
E = Y(2, 1:N);
Ia = Y(3, 1:N);
Iq = Y(4, 1:N);
R = Y(5, 1:N);
D = Y(6, 1:N);
P = Y(7, 1:N);
end