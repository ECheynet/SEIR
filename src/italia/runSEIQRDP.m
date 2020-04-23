clearvars; close all; clc;
addpath('../common/SEIQRDP/');
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

%% setup model
Npop = 60.48e6; % population
Recovered = tableCOVIDItaly_Tot.recovered';
Deaths = tableCOVIDItaly_Tot.dead';
Confirmed = tableCOVIDItaly_Tot.totPositive';
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));

% To simulate the cases after fitting
dt = 1/24; % time step (each hour)
time1 = datetime(time(1)):dt:datetime(datestr(floor(now)+datenum(300)));
N = numel(time1);
t = [0:N-1].*dt;

%% fit
% initial conditions
E0 = 1e-3 * Npop; % exposed
I0 = 1e-2 * Npop; % asymptomatic
Q0 = Confirmed(1);
R0 = Recovered (1);
D0 = Deaths (1);

% initial guess
alpha_guess = 0; % protection rate
beta_guess = 0; % Infection rate
LT_guess = 17; % latent time in days
Q_guess = 0; % rate at which infectious people enter in quarantine
lambda_guess = [0.02, 3.0]; % recovery rate
kappa_guess = [0.02, 0.05]; % death rate
guess = [alpha_guess, beta_guess, 1/LT_guess, Q_guess, lambda_guess, kappa_guess];

% do the fit
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] = fit(Confirmed,Recovered,Deaths,Npop,E0,I0,time,guess,'Display','off');

%% apply model with fitted parameters
[S,E,I,Q,R,D,P] = SEIQRDP(alpha1, beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t);

% errors
x = Q;
[rmseConfirmed, nrmseConfirmed] = mof(Confirmed, x(1:1/dt:length(x)));

x = R;
[rmseRecovered, nrmseRecovered] = mof(Recovered, x(1:1/dt:length(x)));

x = D;
[rmseDeaths, nrmseDeaths] = mof(Deaths, x(1:1/dt:length(x)));

%% plot
figure

semilogy(time1,S/1e1,'c',time1,E/1e1,'r',time1,I/1e1,'b',time1,R/1e1,'k');  % model
hold on

set(gca,'ColorOrderIndex',1);
semilogy(time,Confirmed,'ro',time,Recovered,'bo',time,Deaths,'ko');  % real data

% labels
ylabel('Number of cases')
xlabel('time (days)')
title('Italy');
leg = {'','Confirmed','Recovered','Dead'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')

% prettify
grid on
grid minor
axis tight
set(gca,'yscale','lin')