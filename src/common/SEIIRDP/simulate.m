function [Y] = simulate(alpha, beta, gamma, delta, lambda0, kappa0, tau0, Y, Npop, t)
    N = numel(t);
    dt = median(diff(t));
    iteration = @(Y,A,F) A*Y + F;
    
    % lambda = lambda0(1) * (1 - exp(-lambda0(2).*t));  % todo
    % kappa = kappa0(1) * exp(-kappa0(2).*t);
    kappa = arrayfun(@(x) kappa0(1) * exp(-kappa0(2) * x), t);
    % tau = tau0(1) * (1 - exp(-tau0(2).*t));
        
    for i=1:N-1
        % interactions = [alpha, gamma, delta, lambda(i), kappa(i), tau(i)];  % pack coefficients
        interactions = [alpha, gamma, delta, lambda0(1), kappa(i), tau0(1)];  % pack coefficients
        A = getModelMatrix(interactions);
        SI = Y(1, i) * Y(3, i);  % S * Ia
        F = zeros(7, 1);  % vector just for SI
        F(1:2, 1) = [-beta / Npop; beta / Npop].*SI;
        Y(:, i + 1) = RK4(iteration, Y(:, i), A, F, dt);
    end
end

