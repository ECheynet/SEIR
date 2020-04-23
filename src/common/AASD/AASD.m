function [S,E,I,Q,R,D] = AASD()
% Simulate the time-histories of an epidemic outbreak using AASD (generalized SEIR age stratified) model

%% Initial conditions
N = numel(t);
Y = zeros(6,N);
Y(1,1) = Npop-Q0-E0-R0-D0-I0;  % susceptible
Y(2,1) = E0;  % exposed
Y(3,1) = I0;  % infectious
Y(4,1) = Q0;  % quarantined
Y(5,1) = R0;  % recovered
Y(6,1) = D0;  % dead

%% Matrix version of gSEIR model
[Y] = simulate(beta, gamma, lambda0, kappa0, delta, Y, Npop,t, N);

S = Y(1,1:N);
E = Y(2,1:N);
I = Y(3,1:N);
Q = Y(4,1:N);
R = Y(5,1:N);
D = Y(6,1:N);
end