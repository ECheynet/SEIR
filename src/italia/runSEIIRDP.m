clearvars; close all; clc;
addpath('../common/SEIIRDP/');
addpath('../common/stats/');
addpath('../common/math/');

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
alpha_guess = 0; % protection rate
beta_guess = 0; % S -> E (by coming in contact with asymp)
gamma_guess = 0; % (inverse of latent time in days) rate at which exposed can carry the virus
delta_guess = 0; % asymp -> test positive
lambda_guess = [0.01, 10.0]; % recovery rate (when being symptomatic)
kappa_guess = [1, 0.05]; % death rate (when being symptomatic)
tau_guess = 0;  % asym -> recover
guess = [alpha_guess, beta_guess, gamma_guess, delta_guess, lambda_guess, kappa_guess, tau_guess];

% do the fit
[alpha_fit, beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit] = fit(TotPositive, Recovered, Deaths, Npop, E0, Ia0, time, guess, 'Display', 'off');

%% apply model with fitted parameters
[S, E, Ia, Iq, R, D, P] = model(alpha_fit, beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, Npop, E0, Ia0, Iq0, R0, D0, t);

% errors
x = Iq; % simulated number of total positive
[rmseConfirmed, nrmseConfirmed] = mof(TotPositive, x(1:1/dt:length(x)));

x = R;  % recovered
[rmseRecovered, nrmseRecovered] = mof(Recovered, x(1:1/dt:length(x)));

x = D;  % dead
[rmseDeaths, nrmseDeaths] = mof(Deaths, x(1:1/dt:length(x)));

%% plot
figure
plotter = @plot;

plotter(time1, Iq + R + D); hold on % simulation
plotter(time1, Iq); hold on
plotter(time1, R); hold on
plotter(time1, D); hold on
plotter(time1, E * 1e-1); hold on
plotter(time1, S * 1e-2); hold on
plotter(time1, P * 1e-2); hold on

set(gca,'ColorOrderIndex',1);
plotter(time, TotCases, 'o'); hold on % real data
plotter(time, TotPositive, 'o'); hold on
plotter(time, Recovered, 'o'); hold on
plotter(time, Deaths, 'o');

% labels
ylabel('#')
xlabel('time (days)')
title('Italy');
leg = {'total = positives + recovered + dead', 'positives = quarantined + hospitalized', 'recovered', 'dead', 'exposed * 1e-1', 'susceptible * 1e-2', 'not susceptible * 1e-2'};
legend(leg{:}, 'location', 'southoutside')
set(gcf, 'color', 'w')

% prettify
grid on
grid minor
axis tight
set(gca, 'yscale', 'lin')

%% results summary
latent_period = 1 / gamma_fit
totNRMSE = nrmseConfirmed + nrmseDeaths + nrmseRecovered