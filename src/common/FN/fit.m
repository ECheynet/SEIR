function paramsFit = fit(TotPositive, Recovered, Deaths, Npop, E0, Iu0, time, paramsGuess, varargin)
% Estimates the parameters used in the SEIQRDP function, used to model the time-evolution of an epidemic outbreak.

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX', 1e-5);  %  option for optimset
p.addOptional('tolFun', 1e-2);  %  option for optimset
p.addOptional('Display', 'iter'); % Display option for optimset
p.addOptional('dt', 0.1); % time step for the fitting
p.parse(varargin{:});
tolX = p.Results.tolX;
tolFun = p.Results.tolFun;
Display = p.Results.Display;
dt = p.Results.dt;

%% Options for lsqcurvfit
options = optimset('TolX', tolX, 'TolFun', tolFun, 'MaxFunEvals', 1e6, 'Display', Display);

%% Fitting the data
% Write the target input into a matrix
if size(time,1)>size(time,2) && size(time,2)==1, time = time'; end
fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal 
t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges

% lsqcurvefit settings
guess = paramsGuess.getAsVector();
lowerBounds = paramsGuess.getLowerBounds();
upperBounds = paramsGuess.getUpperBounds();
f = @(params, t) optim(params, t);
[Coeff,~,~,~,~,~,~] = lsqcurvefit(f, guess, tTarget(:)', [TotPositive; Recovered; Deaths], lowerBounds, upperBounds, options);

%% Write the fit coeff in the outputs
paramsFit = [ ...
    Coeff(1), ...
    Coeff(2), ...
    Coeff(3), ...
    Coeff(4), ...
    Coeff(5), ...
    Coeff(6), ...
    Coeff(7), ...
    Coeff(8) ...
];

    %% optimization function
    function [output] = optim(params, t0)
        %% Initial conditions
        N = numel(t);
        Y = zeros(7, N);

        Y(1,1) = Npop - TotPositive(1) - Recovered(1) - Deaths(1) - E0 - Iu0;
        Y(2,1) = E0;
        Y(3,1) = Iu0;
        Y(4,1) = TotPositive(1);
        Y(5,1) = Recovered(1);
        Y(6,1) = Deaths(1);
        Y(7,1) = 0; % P state

        %% ODE solution
        [Y] = simulate(params, Y, Npop, t);

        Q1 = Y(4,1:N); % confirmed
        R1 = Y(5,1:N); % recovered
        D1 = Y(6,1:N); % dead

        Q1 = interp1(t, Q1, t0);
        R1 = interp1(t, R1, t0);
        D1 = interp1(t, D1, t0);

        output = [Q1; R1; D1];
    end
end