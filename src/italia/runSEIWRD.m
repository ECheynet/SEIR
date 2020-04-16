clearvars; close all; clc;
addpath('../common/SEIQRD/');
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
Recovered = tableCOVIDItaly_Tot.Recovered';
Deaths = tableCOVIDItaly_Tot.Deaths';
Confirmed = tableCOVIDItaly_Tot.Confirmed';
Quarantined = tableCOVIDItaly_Tot.Quarantined';
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));

% if the number of confirmed Confirmed cases is small, it is difficult to know whether the quarantine has been rigorously applied or not. In addition, this suggests that the number of infectious is much larger than the number of confirmed cases
minNum = round(0.1*max(Confirmed));
Recovered(Confirmed<=minNum) = [];
Deaths(Confirmed<=minNum) = [];
time(Confirmed<=minNum) = [];
Quarantined(Confirmed<=minNum) = [];
Confirmed(Confirmed<=minNum) = [];

% To simulate the cases after fitting
dt = 1/24; % time step (each hour)
time1 = datetime(time(1)):dt:datetime(datestr(floor(now)+datenum(7)));
N = numel(time1);
t = [0:N-1].*dt;

%% fit
% initial conditions
E0 = Confirmed(1); % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 = Confirmed(1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = Quarantined(1);
R0 = Recovered(1);
D0 = Deaths(1);

% initial guess
beta_guess = 3.0; % Infection rate
LT_guess = 22; % latent time in days
Q_guess = 0.24; % rate at which infectious people enter in quarantine
lambda_guess = [0.02, 3.0]; % recovery rate
kappa_guess = [0.02, 0.05]; % death rate
guess = [beta_guess, 1/LT_guess, Q_guess, lambda_guess, kappa_guess];

% do the fit
[beta1,gamma1,delta1,Lambda1,Kappa1] = fit(Quarantined,Recovered,Deaths,Npop,E0,I0,time,guess,'Display','off');

%% apply model with fitted parameters
[S,E,I,Q,R,D] = SEIQRD(beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t);

% errors
x = Q+R+D;
[rmseConfirmed, nrmseConfirmed] = mof(Confirmed, x(1:1/dt:length(x)));

x = R;
[rmseRecovered, nrmseRecovered] = mof(Recovered, x(1:1/dt:length(x)));

x = D;
[rmseDeaths, nrmseDeaths] = mof(Deaths, x(1:1/dt:length(x)));

%% plot
figure

semilogy(time1,Q+R+D,'c',time1,Q,'r',time1,R,'b',time1,D,'k');  % model
hold on

set(gca,'ColorOrderIndex',1);
semilogy(time,Confirmed,'co',time,Quarantined,'ro',time,Recovered,'bo',time,Deaths,'ko');  % real data

% labels
ylabel('Number of cases')
xlabel('time (days)')
title('Italy');
leg = {'Total','Quarantined (confirmed infectious)','Recovered','Dead'};
legend(leg{:},'location','southoutside')
set(gcf,'color','w')

% prettify
grid on
grid minor
axis tight
set(gca,'yscale','lin')