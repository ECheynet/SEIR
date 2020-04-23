function [beta_fit, gamma_fit, delta_fit, lambda_fit, kappa_fit, tau_fit, rho_fit] = fit(TotPositive, Recovered, Deaths, Npop, E0, Ia0, time, guess, varargin)
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
upperBounds = 10 * ones(1, numel(guess));
[Coeff,~,~,~,~,~,~] = lsqcurvefit(@(para,t) SEIIRD_for_fitting(para, t), guess, tTarget(:)', [TotPositive; Recovered; Deaths], lowerBounds, upperBounds, options);

%% Write the fitted coeff in the outputs
beta_fit = abs(Coeff(1));
gamma_fit = abs(Coeff(2));
delta_fit = abs(Coeff(3));
lambda_fit = abs(Coeff(4:5));
kappa_fit = abs(Coeff(6:7));
tau_fit = abs(Coeff(7:8));
rho_fit = abs(Coeff(8:9));

%% nested functions
    function [output] = SEIIRD_for_fitting(para, t0)
        beta = abs(para(1));
        gamma = abs(para(2));
        delta = abs(para(3));
        lambda0 = abs(para(4:5));
        kappa0 = abs(para(6:7));
        tau0 = abs(para(7:8));
        rho0 = abs(para(8:9));
        
        %% Initial conditions
        N = numel(t);
        Y = zeros(6, N);
        
        Y(1,1) = Npop - TotPositive(1) - Recovered(1) - Deaths(1) - E0 - Ia0;
        Y(2,1) = E0;
        Y(3,1) = Ia0;
        Y(4,1) = TotPositive(1);
        Y(5,1) = Recovered(1);
        Y(6,1) = Deaths(1);
        
        if round(sum(Y(:,1)) - Npop) ~= 0
            error('the sum must be zero because the total population (including the deads) is assumed constant');
        end
        
        %% ODE solution
        [Y] = simulate(beta, gamma, delta, lambda0, kappa0, tau0, rho0, Y, Npop, t);
        
        Q1 = Y(4,1:N);
        R1 = Y(5,1:N);
        D1 = Y(6,1:N);
        
        Q1 = interp1(t, Q1, t0);
        R1 = interp1(t, R1, t0);
        D1 = interp1(t, D1, t0);
        
        output = [Q1; R1; D1];
    end
end