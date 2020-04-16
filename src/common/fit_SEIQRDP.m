function [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,varargout] = fit_SEIQRDP(Q,R,D,Npop,E0,I0,time,guess,varargin)
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,varargout] = 
% fit_SEIQRDP(Q,R,D,Npop,E0,I0,time,guess,varargin) estimates the 
% parameters used in the SEIQRDP function, used to model the time-evolution
% of an epidemic outbreak.
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
%   alpha: scalar [1x1]: fitted protection rate
%   beta: scalar [1x1]: fitted  infection rate
%   gamma: scalar [1x1]: fitted  Inverse of the average latent time
%   delta: scalar [1x1]: fitted  inverse of the average quarantine time
%   lambda: scalar [1x1]: fitted  cure rate
%   kappa: scalar [1x1]: fitted  mortality rate
%   optional:
%       - residual
%       - Jcobian
%       - The function @SEIQRDP_for_fitting
% 
% Author: E. Cheynet - UiB - last modified 24-03-2020
% 
% see also SEIQRDP.m

%%

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('tolX',1e-4);  %  option for optimset
p.addOptional('tolFun',1e-4);  %  option for optimset
p.addOptional('Display','iter'); % Display option for optimset
p.addOptional('dt',0.1); % time step for the fitting
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if size(time,1)>size(time,2) && size(time,2)==1,    time = time';end
if size(time,1)>1 && size(time,2)>1,  error('Time should be a vector');end

fs = 1./dt;
tTarget = round(datenum(time-time(1))*fs)/fs; % Number of days with one decimal 

t = tTarget(1):dt:tTarget(end); % oversample to ensure that the algorithm converges



modelFun1 = @SEIQRDP_for_fitting; % transform a nested function into anonymous function

% call Lsqcurvefit
[Coeff,~,residual,~,~,~,jacobian] = lsqcurvefit(@(para,t) modelFun1(para,t),...
    guess,tTarget(:)',input,zeros(1,numel(guess)),[1 3 1 1 2 3 2 2],options);


if nargout ==7
    varargout{1} = residual;
elseif nargout==8
    varargout{1} = residual;
    varargout{2} = jacobian;
elseif nargout==9
    varargout{1} = residual;
    varargout{2} = jacobian;
    varargout{3} = modelFun1;
elseif nargout>9
    error('Too many output specified')
end


%% Write the fitted coeff in the outputs

alpha1 = abs(Coeff(1));
beta1 = abs(Coeff(2));
gamma1 = abs(Coeff(3));
delta1 = abs(Coeff(4));
Lambda1 = abs(Coeff(5:6));
Kappa1 = abs(Coeff(7:8));

% if isempty(R)
%     Lambda1(2)=0;
% end










%% nested functions

    function [output] = SEIQRDP_for_fitting(para,t0)

        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5:6));
        kappa0 = abs(para(7:8));

        
        %% Initial conditions
        N = numel(t);
        Y = zeros(7,N);
        Y(2,1) = E0;
        Y(3,1) = I0;
        Y(4,1) = Q(1);
        if ~isempty(R),
            Y(5,1) = R(1);
            Y(1,1) = Npop-Q(1)-R(1)-D(1)-E0-I0;
        else
            Y(1,1) = Npop-Q(1)-D(1)-E0-I0;
        end
        Y(6,1) = D(1);
        
        if round(sum(Y(:,1))-Npop)~=0
            error('the sum must be zero because the total population (including the deads) is assumed constant');
        end
        %%
        modelFun = @(Y,A,F) A*Y + F;
        
         
         
%          if ~isempty(R)
                  lambda = lambda0(1)*(1-exp(-lambda0(2).*t)); % I use these functions for illustrative purpose only
                  kappa = kappa0(1)*exp(-kappa0(2).*t); 
%          else
%              lambda = lambda0(1).*ones(1,N); % I use these functions for illustrative purpose only
%              kappa = kappa0(1)*exp(-kappa0(2).*t); 
%          end
                  
        % ODE reYution
        for ii=1:N-1
            A = getA(alpha,gamma,delta,lambda(ii),kappa(ii));
            SI = Y(1,ii)*Y(3,ii);
            F = zeros(7,1);
            F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
            Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
        end
        
%         I1 = Y(3,1:N);
        Q1 = Y(4,1:N);
        R1 = Y(5,1:N);
        D1 = Y(6,1:N);
        
        Q1 = interp1(t,Q1,t0);
        R1 = interp1(t,R1,t0);
        D1 = interp1(t,D1,t0);
        if ~isempty(R),
            output = [Q1;R1;D1];
        else
            output = [Q1+R1;D1];
        end
        
    end
    function [A] = getA(alpha,gamma,delta,lambda,kappa)
        
        A = zeros(7);
        % S
        A(1,1) = -alpha;
        % E
        A(2,2) = -gamma;
        % I
        A(3,2:3) = [gamma,-delta];
        % Q
        A(4,3:4) = [delta,-kappa-lambda];
        % R
        A(5,4) = lambda;
        % D
        A(6,4) = kappa;
        % P
        A(7,1) = alpha;
        
    end
    function [Y] = RK4(Fun,Y,A,F,dt)
        
        % Runge-Kutta of order 4
        k_1 = Fun(Y,A,F);
        k_2 = Fun(Y+0.5*dt*k_1,A,F);
        k_3 = Fun(Y+0.5*dt*k_2,A,F);
        k_4 = Fun(Y+k_3*dt,A,F);
        % output
        Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt;
    end

end

