clearvars; close all; clc;
addpath('../common/FN/');
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
undetectedDeaths = 0.36;  % see Gabanelli in Corriere
UndetectedDeaths = Deaths * undetectedDeaths;
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));

% To simulate the cases after fitting
dt = 1/24; % time step (each hour)
daysToPredict = 2 * 30;  % 1 month
time1 = datetime(time(1)) : dt : datetime(datestr(floor(now) + datenum(daysToPredict)));
N = numel(time1);
t = [0:N - 1].*dt;

%% fit
% initial conditions
E0 = 1e-3 * Npop; % starting exposed = 0.1% population
Iu0 = 0.3 * Npop; % asymptomatic = 30% population (see Vo' study: 43 %)
Iq0 = TotPositive(1);
R0 = Recovered(1);
D0 = Deaths(1);
C0 = Deaths(1) * undetectedDeaths;

% initial guess
alpha_guess = 0; % protection rate
beta_guess = 0; % S -> E (by coming in contact with asymp)
gamma_guess = 1/17; % (inverse of latent time in days) rate at which exposed can carry the virus
delta_guess = 0; % asymp -> test positive
lambda_guess = [0.01, 1]; % recovery rate (when being symptomatic)
kappa_guess = [1, 0.01]; % death rate (when being symptomatic)
tau_guess = [0.01, 1];  % asym -> recover
rho_guess = [1, 0.01]; % death rate (when being asymptomatic)
guess = [alpha_guess, beta_guess, gamma_guess, delta_guess, lambda_guess, kappa_guess, tau_guess, rho_guess];

% do the fit
[alpha_fit, beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, rho_fit] = fit(TotPositive, Recovered, Deaths, UndetectedDeaths, Npop, E0, Iu0, time, guess, 'Display', 'off');

%% apply model with fitted parameters
[S, E, Iu, Iq, R, D, P, C] = model(alpha_fit, beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, rho_fit, Npop, E0, Iu0, Iq0, R0, D0, C0, t);

% errors
x = Iq; % simulated number of total positive
[rmseConfirmed, nrmseConfirmed] = mof(TotPositive, x(1:1/dt:length(x)));

x = R;  % recovered
[rmseRecovered, nrmseRecovered] = mof(Recovered, x(1:1/dt:length(x)));

x = D + C;  % dead
[rmseDeaths, nrmseDeaths] = mof(Deaths + UndetectedDeaths, x(1:1/dt:length(x)));

%% plot
figure
plotter = @plot;

yyaxis left  % use left y-axis
p1 = plotter(time1, Iq + R + D, '-b'); hold on % simulation
p2 = plotter(time1, Iq, '-r'); hold on
p3 = plotter(time1, R, '-g'); hold on
p4 = plotter(time1, D + C, '-k'); hold on

p5 = plotter(time, TotCases, 'xb'); hold on % real data
p6 = plotter(time, TotPositive, 'xr'); hold on
p7 = plotter(time, Recovered, 'xg'); hold on
p8 = plotter(time, Deaths + UndetectedDeaths, 'xk');

yyaxis right  % use right y-axis for large numbers
p9 = plotter(time1, E, '-c'); hold on
p10 = plotter(time1, S, '--k'); hold on
p11 = plotter(time1, P, ':k'); hold on
p12 = plotter(time1, Iu, '-.k'); hold on

% labels
ylabel('#')
xlabel('time (days)')
title('Italy');
leg = {'total cases = positives + recovered + dead', ...
    'total positives = quarantined + hospitalized', ...
    'total recovered', ...
    'total dead', ...
    'real total cases', ...
    'real total positives', ...
    'real total recovered', ...
    'real total dead (undetected too)', ...
    'exposed', ...
    'susceptible', ...
    'not susceptible', ...
    'asymptomatic'};
legend([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12], leg{:});
set(gcf, 'color', 'w')

% prettify
grid on
grid minor
axis tight

%% results summary
latent_period = 1 / gamma_fit
totNRMSE = nrmseConfirmed + nrmseDeaths + nrmseRecovered