clearvars; close all; clc;
addpath('../common/SEIIRD/');
addpath('../common/stats/');

%% get data
tableCOVIDItaly = getData();
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));
fprintf(['Most recent update: ',datestr(time(end)),'\n'])

% parse data
% merge regional data for each day
tableCOVIDItaly_Tot = varfun(@sum,tableCOVIDItaly, 'InputVariables',tableCOVIDItaly.Properties.VariableNames(7:end), 'GroupingVariables','Date');

% remove the 'GroupCount' variable, should total to the number of Italian regions (19 + 2 autonomous provinces)
tableCOVIDItaly_Tot = removevars(tableCOVIDItaly_Tot,'GroupCount');

% rename the accumulated variables with the original variable names
tableCOVIDItaly_Tot.Properties.VariableNames=[tableCOVIDItaly.Properties.VariableNames(1),tableCOVIDItaly.Properties.VariableNames(7:end)];

Recovered = tableCOVIDItaly_Tot.recovered';
Deaths = tableCOVIDItaly_Tot.dead';
% Hospitalized = tableCOVIDItaly_Tot.hospitalized';
% Quarantined = tableCOVIDItaly_Tot.quarantined';
TotPositive = tableCOVIDItaly_Tot.totPositive'; % = #quarantined + #hospitalized
TotCases = tableCOVIDItaly_Tot.totCases'; % = #totPositive + #recovered + #dead

%% setup model
Npop = 60.48e6; % population
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));

% To simulate the cases after fitting
dt = 1/24; % time step (each hour)
daysToPredict = 14;
time1 = datetime(time(1)) : dt : datetime(datestr(floor(now) + datenum(daysToPredict)));
N = numel(time1);
t = [0:N - 1].*dt;

%% fit
% initial conditions
E0 = 1e-3 * Npop; % exposed
Ia0 = 1e-2 * Npop; % asymptomatic
Iq0 = TotPositive(1);
R0 = Recovered(1);
D0 = Deaths(1);

% initial guess
beta_guess = 0; % Infection rate
gamma_guess = 1 / 17; % (latent time in days) rate at which exposed can carry the virus
delta_guess = 0; % rate at which exposed go confirmed (and quarantined)
lambda_guess = zeros(1, 2); % recovery rate (when being asymptomatic)
kappa_guess = zeros(1, 2); % death rate (when being asymptomatic)
tau_guess = zeros(1, 2); % recovery rate (when being confirmed)
rho_guess = zeros(1, 2); % death rate (when being confirmed)
guess = [beta_guess, gamma_guess, delta_guess, lambda_guess, kappa_guess, tau_guess, rho_guess];

% do the fit
[beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, rho_fit] = fit(TotPositive, Recovered, Deaths, Npop, E0, Ia0, time, guess, 'Display', 'off');

%% apply model with fitted parameters
[S, E, Ia, Iq, R, D] = SEIIRD(beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, rho_fit, Npop, E0, Ia0, Iq0, R0, D0, t);

% errors
x = Iq; % simulated number of total positive
[rmseConfirmed, nrmseConfirmed] = mof(TotPositive, x(1:1/dt:length(x)));

x = R;  % recovered
[rmseRecovered, nrmseRecovered] = mof(Recovered, x(1:1/dt:length(x)));

x = D;  % dead
[rmseDeaths, nrmseDeaths] = mof(Deaths, x(1:1/dt:length(x)));

%% plot
figure

semilogy(time1, E, 'v', time1, Ia, 'blue', time1, Iq + R + D, 'cyan', time1, Iq, 'red', time1, R, 'green', time1, D, 'black');  % simulation
hold on

set(gca,'ColorOrderIndex',1);
semilogy(time, TotCases, 'co', time, TotPositive, 'ro', time, Recovered, 'go', time, Deaths, 'ko');  % real data

% labels
ylabel('#')
xlabel('time (days)')
title('Italy');
leg = {'exposed', 'asymptomatic', 'total = positives + recovered + dead', 'positives = quarantined + hospitalized', 'recovered', 'dead'};
legend(leg{:}, 'location', 'southoutside')
set(gcf, 'color', 'w')

% prettify
grid on
grid minor
axis tight
set(gca, 'yscale', 'lin')