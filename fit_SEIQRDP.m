function [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] = fit_SEIQRDP(I,R,D,Npop,time,guess)
% [alpha1,beta1,gamma1,delta1,Lambda1,Kappa1] =
% fit_SEIQRDP(I,R,D,Npop,time,guess) estimates the parameters used in the
% SEIQRDP function, used to model the time-evolution of an epidemic outbreak.
% 
% Input
% 
%   I: vector [1xN] of the target time-histories of the infectious cases
%   R: vector [1xN] of the target time-histories of the recovered cases
%   D: vector [1xN] of the target time-histories of the dead cases
%   Npop: scalar: Total population of the sample
%   time: vector [1xN] of time (datetime)
%   guess: first vector [1x6] guess for the fit
% 
% Output
% 
%   alpha: scalar [1x1]: fitted protection rate
%   beta: scalar [1x1]: fitted  infection rate
%   gamma: scalar [1x1]: fitted  Inverse of the average latent time
%   delta: scalar [1x1]: fitted  inverse of the average quarantine time
%   lambda: scalar [1x1]: fitted  cure rate
%   kappa: scalar [1x1]: fitted  mortality rate
% 
% Author: E. Cheynet - UiB - last modified 16-03-2020
% 
% see also SEIQRDP.m

%%

options=optimset('TolX',1e-4,'TolFun',1e-4,'Display','iter');
input = [I;R;D];
tTarget = datenum(time-time(1)); % Number of days
t = tTarget(1):0.1:tTarget(end);
dt = median(diff(t));


modelFun1 = @SEIQRDP_for_fitting; % transform a nested function into anonymous function
Coeff = lsqcurvefit(@(para,t) modelFun1(para,t),guess,tTarget(:)',input,[0,0,0,0,0,0],[3,3,2,2,3,3],options);

alpha1 = abs(Coeff(1));
beta1 = abs(Coeff(2));
gamma1 = abs(Coeff(3));
delta1 = abs(Coeff(4));
Lambda1 = abs(Coeff(5));
Kappa1 = abs(Coeff(6));


    function [output] = SEIQRDP_for_fitting(para,t0)

        alpha = abs(para(1));
        beta = abs(para(2));
        gamma = abs(para(3));
        delta = abs(para(4));
        lambda0 = abs(para(5));
        kappa0 = abs(para(6));

        
        %% Initial conditions
        N = numel(t);
        Y = zeros(7,N);
        Y(1,1) = Npop-I(1)-R(1)-D(1);
        Y(3,1) = I(1);
        Y(5,1) = R(1);
        Y(6,1) = D(1);
        if round(sum(Y(:,1))-Npop)~=0
            error('the sum must be zero because the total population (including the deads) is assumed constant');
        end
        %%
        modelFun = @(Y,A,F) A*Y + F;
        dt = median(diff(t));
        % ODE reYution
        for ii=1:N-1
            
            lambda = 0.03*(1-exp(-lambda0.*t(ii))); % I use these functions for illustrative purpose only
            kappa = 0.06*exp(-kappa0.*t(ii)); % I use these functions for illustrative purpose only
            A = getA(alpha,gamma,delta,lambda,kappa);
            SI = Y(1,ii)*Y(3,ii);
            F = zeros(7,1);
            F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
            Y(:,ii+1) = RK4(modelFun,Y(:,ii),A,F,dt);
        end
        
        
        I0 = Y(3,1:N);
        R0 = Y(5,1:N);
        D0 = Y(6,1:N);
        
        I0 = interp1(t,I0,t0);
        R0 = interp1(t,R0,t0);
        D0 = interp1(t,D0,t0);
        
        output = [I0;R0;D0];
%         output = I0;
        
        
        
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

