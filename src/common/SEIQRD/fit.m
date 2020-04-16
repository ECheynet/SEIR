function [beta1,gamma1,delta1,Lambda1,Kappa1] = fit(Q,R,D,Npop,E0,I0,time,guess,varargin)
% Estimates the parameters used in the SEIQRDP function, used to model the time-evolution of an epidemic outbreak.
% 
% Input
% 
%   I: vector [1xN] of the target time-histories of the infectious cases
%   R: vector [1xN] of the target time-histories of the recovered cases
%   D: vector [1xN] of the target time-histories of the dead cases
%   Npop: scalar: Total population of the sample
%   E0: scalar [1x1]: Initial number of exposed cases
%   I0: scalar [1x1]: Initial number of infectious cases
%   time: vector [1xN] of time (datetime)
%   guess: first vector [1x6] guess for the fit
%   optionals
%       -tolFun: tolerance  option for optimset
%       -tolX: tolerance  option for optimset
%       -Display: Display option for optimset
%       -dt: time step for the fitting function
% 
% Output
% 
%   beta: scalar [1x1]: fitted  infection rate
%   gamma: scalar [1x1]: fitted  Inverse of the average latent time
%   delta: scalar [1x1]: fitted  inverse of the average quarantine time
%   lambda: scalar [1x1]: fitted  cure rate
%   kappa: scalar [1x1]: fitted  mortality rate

%%

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  %  option for optimset
p.addOptional('tolFun',1e-4);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting
p.parse(varargin{:});
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
Display  = p.Results.Display ;
dt  = p.Results.dt ;

%% Options for lsqcurvfit
options=optimset('TolX',tolX,'TolFun',tolFun,'MaxFunEvals',800,'Display',Display);

%% Fitting the data
% Write the target input into a matrix
Q(Q<0)=0; % negative values are not possible
R(R<0)=0; % negative values are not possible
D(D<0)=0; % negative values are not possible

if isempty(R)
    warning(' No data available for "Recovered" ')
    input = [Q;D];
else
    input = [Q;R;D];
end

if size(time,1)>size(time,2) && size(time,2)==1, time = time'; end
if size(time,1)>1 && size(time,2)>1, error('Time should be a vector'); end

fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal 

t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges

modelFun1 = @SEIQRD_for_fitting; % transform a nested function into anonymous function

% call Lsqcurvefit
[Coeff,~,~,~,~,~,~] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess,tTarget(:)',input,zeros(1,numel(guess)),[1 3 1 1 2 3 2 2],options);

%% Write the fitted coeff in the outputs
beta1 = abs(Coeff(1));
gamma1 = abs(Coeff(2));
delta1 = abs(Coeff(3));
Lambda1 = abs(Coeff(4:5));
Kappa1 = abs(Coeff(6:7));

%% nested functions
    function [output] = SEIQRD_for_fitting(para,t0)
        beta = abs(para(1));
        gamma = abs(para(2));
        delta = abs(para(3));
        lambda0 = abs(para(4:5));
        kappa0 = abs(para(6:7));
        
        %% Initial conditions
        N = numel(t);
        Y = zeros(6,N);
        Y(2,1) = E0;
        Y(3,1) = I0;
        Y(4,1) = Q(1);
        if ~isempty(R)
            Y(5,1) = R(1);
            Y(1,1) = Npop-Q(1)-R(1)-D(1)-E0-I0;
        else
            Y(1,1) = Npop-Q(1)-D(1)-E0-I0;
        end
        Y(6,1) = D(1);
        
        if round(sum(Y(:,1))-Npop)~=0
            error('the sum must be zero because the total population (including the deads) is assumed constant');
        end
        
        %% ODE solution
        [Y] = simulate(beta, gamma, lambda0, kappa0, delta, Y, Npop, t, N);
        
        Q1 = Y(4,1:N);
        R1 = Y(5,1:N);
        D1 = Y(6,1:N);
        
        Q1 = interp1(t,Q1,t0);
        R1 = interp1(t,R1,t0);
        D1 = interp1(t,D1,t0);
        if ~isempty(R)
            output = [Q1;R1;D1];
        else
            output = [Q1+R1;D1];
        end
    end
end

