function [Y] = simulate(beta, gamma, delta, lambda0, kappa0, tau0, rho0, Y, Npop, t)
    N = numel(t);
    dt = median(diff(t));
    iteration = @(Y,A,F) A*Y + F;
    
    lambda = lambda0(1) * (1 - exp(-lambda0(2).*t));  % todo (illustrative purpose only)
    kappa = kappa0(1) * exp(-kappa0(2).*t);
    tau =  tau0(1) * (1 - exp(-tau0(2).*t));
    rho = rho0(1) * exp(-rho0(2).*t);
        
    for i=1:N-1
        A = getIterationMatrix(gamma, delta, lambda(i), kappa(i), tau(i), rho(i));
        SI = Y(1, i) * Y(3, i);  % S * Ia
        F = zeros(6, 1);  % vector just for SI
        F(1:2, 1) = [-beta / Npop; beta / Npop].*SI;
        Y(:, i + 1) = RK(iteration, Y(:, i), A, F, dt);
    end
end

