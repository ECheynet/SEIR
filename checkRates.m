function checkRates(time,Q,R,D,kappaFun,lambdaFun,kappa,lambda)
% function checkRates(time,Q,D,R,kappaFun,lambdaFun) compares the fitted
% and calcualted death and recovered ratios. The idea is to check whether
% the approximation of these ratios is appropriate
%
% Inputs
%   time: datetime: [1xN]: time array
%   Q: double [1xN]: Time histories of the quarantined/active cases :
%   D: double [1xN]: Time histories of the deceased cases :
%   R: double [1xN]: Time histories of the recovered cases :
%   kappaFun: anonymous function approximating the death rate
%   lambdaFun: anonymous function approximating the recovery rate
%
% Outputs:
% None
%
% Author: E. Cheynet - UiB - last modified 07-05-2020
%
% see also SEIQRDP.m fit_SEIQRDP.m

%% Compute the rate of deceased and recovered cases

Q = Q(:);
R = R(:);
D = D(:);
time = time(:);

rateD = (diff(D)./diff(datenum(time-time(1))))./Q(2:end);
rateD(abs(rateD)>3) = nan; % remove bovious outliers

if ~isempty(R)
    rateR = (diff(R)./diff(datenum(time-time(1))))./Q(2:end);
    rateR(abs(rateR)>3) = nan;
end
%% Define the time
x = datenum(time(2:end)-time(1));
x1 = x(1):1/24:x(end);

%% Compare the fitted and effective rates

if ~isempty(R)
    figure;
    subplot(121)
    title('Death rate')
    plot(x,rateD,'k*',x1,kappaFun(kappa,x1),'r')
    xlabel('Time (days)')
    ylabel('Death rate (day^{-1})')
    axis tight
    legend('Measured','Fitted')
    
    
    subplot(122)
    title('Recovery rate')
    
    plot(x,rateR,'b*',x1,lambdaFun(lambda,x1),'r')
    axis tight
    set(gcf,'color','w')
    xlabel('Time (days)')
    ylabel('Recovery rate (day^{-1})')
    legend('Measured','Fitted')
else
    figure;
    plot(x,rateD,'k*',x1,kappaFun(kappa,x1),'r')
    xlabel('Time (days)')
    ylabel('Pseudo-death rate (day^{-1})')
    axis tight
    legend('Measured','Fitted')
end
end

