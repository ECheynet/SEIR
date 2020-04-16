function [S,E,I,Q,R,D] = SEIQRD(beta,gamma,delta,lambda0,kappa0,Npop,E0,I0,Q0,R0,D0,t)
% Simulate the time-histories of an epidemic outbreak using gSEIR (generalized SEIR) model
%
% Input
%   beta: scalar [1x1]: infection rate
%   gamma: scalar [1x1]: inverse of the average latent time
%   delta: scalar [1x1]: inverse of the average quarantine time
%   lambda: scalar [1x1]: cure rate
%   kappa: scalar [1x1]: mortality rate
%   Npop: scalar: Total population of the sample
%   E0: scalar [1x1]: Initial number of exposed cases
%   I0: scalar [1x1]: Initial number of infectious cases
%   Q0: scalar [1x1]: Initial number of quarantined cases
%   R0: scalar [1x1]: Initial number of recovered cases
%   D0: scalar [1x1]: Initial number of dead cases
%   t: vector [1xN] of time (double; it cannot be a datetime)
%
% Output
%   S: vector [1xN] of the target time-histories of the susceptible cases
%   E: vector [1xN] of the target time-histories of the exposed cases
%   I: vector [1xN] of the target time-histories of the infectious cases
%   Q: vector [1xN] of the target time-histories of the quarantinedcases
%   R: vector [1xN] of the target time-histories of the recovered cases
%   D: vector [1xN] of the target time-histories of the dead cases

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