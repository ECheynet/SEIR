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
TotPositive = tableCOVIDItaly_Tot.totPositive'; % = #quarantined + #hospitalized
TotCases = tableCOVIDItaly_Tot.totCases'; % = #totPositive + #recovered + #dead

%% setup model
Npop = 60.48e6; % population
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));

% To simulate the cases after fitting
dt = 1/24; % time step (each hour)
daysToPredict = 30;
time1 = datetime(time(1)) : dt : datetime(datestr(floor(now) + datenum(daysToPredict)));
N = numel(time1);
t = [0:N - 1].*dt;

%% fit
% initial conditions
E0 = 0; % Initial number of exposed cases. Unknown but unlikely to be zero.
Ia0 = 0.1 * TotPositive(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Iq0 = TotPositive(1);
R0 = Recovered(1);
D0 = Deaths(1);

% initial guess
beta_guess = 4; % Infection rate
gamma_guess = 1 / 4; % (latent time in days) rate at which exposed go asymptomatic
delta_guess = 1; % rate at which exposed go confirmed (and quarantined)
lambda_guess = [0.5, 0.05]; % recovery rate (when being asymptomatic)
kappa_guess = [0.1, 0.01]; % death rate (when being asymptomatic)
tau_guess = [0.01, 0.02]; % recovery rate (when being confirmed)
rho_guess = [0.02, 0.02]; % death rate (when being confirmed)
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

semilogy(time1, Iq + R + D, 'c', time1, Iq, 'r', time1, R, 'b', time1, D, 'k');  % simulation
hold on

set(gca,'ColorOrderIndex',1);
semilogy(time, TotCases, 'co', time, TotPositive, 'ro', time, Recovered, 'bo', time, Deaths, 'ko');  % real data

% labels
ylabel('# cases')
xlabel('time (days)')
title('Italy');
leg = {'total = # positives + # recovered + # dead', '# positives = # quarantined + # hospitalized', '# recovered', '# dead'};
legend(leg{:}, 'location', 'southoutside')
set(gcf, 'color', 'w')

% prettify
grid on
grid minor
axis tight
set(gca, 'yscale', 'lin')