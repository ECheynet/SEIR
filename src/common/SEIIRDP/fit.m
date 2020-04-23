function [alpha_fit, beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit] = fit(TotPositive, Recovered, Deaths, Npop, E0, Ia0, time, guess, varargin)
% Estimates the parameters used in the SEIQRDP function, used to model the time-evolution of an epidemic outbreak.

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX', 1e-4);  %  option for optimset
p.addOptional('tolFun', 1e-4);  %  option for optimset
p.addOptional('Display', 'iter'); % Display option for optimset
p.addOptional('dt', 0.1); % time step for the fitting
p.parse(varargin{:});
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display = p.Results.Display ;
dt = p.Results.dt ;

%% Options for lsqcurvfit
options = optimset('TolX', tolX, 'TolFun', tolFun, 'MaxFunEvals', 800, 'Display', Display);

%% Fitting the data
% Write the target input into a matrix
if size(time,1)>size(time,2) && size(time,2)==1, time = time'; end
if size(time,1)>1 && size(time,2)>1, error('Time should be a vector'); end

fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal 
t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges

% call Lsqcurvefit
lowerBounds = zeros(1, numel(guess));
upperBounds = 5 * ones(1, numel(guess));
[Coeff,~,~,~,~,~,~] = lsqcurvefit(@(para,t) optim(para, t), guess, tTarget(:)', [TotPositive; Recovered; Deaths], lowerBounds, upperBounds, options);

%% Write the fitted coeff in the outputs
alpha_fit = abs(Coeff(1));
beta_fit = abs(Coeff(2));
gamma_fit = abs(Coeff(3));
delta_fit = abs(Coeff(4));
lambda_fit = abs(Coeff(5:6));
kappa_fit = abs(Coeff(7:8));
tau_fit = abs(Coeff(9));

%% nested functions
    function [output] = optim(para, t0)
        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5:6));
        kappa0 = abs(para(7:8));
        tau = abs(para(9));
        
        %% Initial conditions
        N = numel(t);
        Y = zeros(7, N);
        
        Y(1,1) = Npop - TotPositive(1) - Recovered(1) - Deaths(1) - E0 - Ia0;
        Y(2,1) = E0;
        Y(3,1) = Ia0;
        Y(4,1) = TotPositive(1);
        Y(5,1) = Recovered(1);
        Y(6,1) = Deaths(1);
        Y(7,1) = 0;
        
        %% ODE solution
        [Y] = simulate(alpha, beta, gamma, delta, lambda0, kappa0, tau, Y, Npop, t);
        
        Q1 = Y(4,1:N);
        R1 = Y(5,1:N);
        D1 = Y(6,1:N);
        
        Q1 = interp1(t, Q1, t0);
        R1 = interp1(t, R1, t0);
        D1 = interp1(t, D1, t0);
        
        output = [Q1; R1; D1];
    end
end