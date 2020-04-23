function [Y] = simulate(beta, gamma, lambda0, kappa0, delta, Y, Npop,t,N)
    dt = median(diff(t));
    iteration = @(Y,A,F) A*Y + F;
    lambda = lambda0(1)*(1-exp(-lambda0(2).*t));  % todo (illustrative purpose only)
    kappa = kappa0(1)*exp(-kappa0(2).*t); 
        
    for ii=1:N-1
        A = getA(gamma, delta, lambda(ii), kappa(ii));
        SI = Y(1,ii)*Y(3,ii);
        F = zeros(6, 1);
        F(1:2,1) = [-beta/Npop;beta/Npop].*SI;
        Y(:,ii+1) = RK4(iteration,Y(:,ii),A,F,dt);
    end
end

