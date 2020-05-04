function [Y] = simulate(paramsFit, Y, Npop, t)
    N = numel(t);
    dt = median(diff(t));
    iteration = @(Y,A,F) A*Y + F;

    lambda0 = paramsFit.lambda;
    lambda = lambda0(1) * (1 - exp(-lambda0(2).*t));  % exponential growth of recovery methods
    
    for i=1:N-1
        paramsAtIteration = paramsFit;
        paramsAtIteration.lambda = lambda(i);  % pack coefficients
        
        A = paramsAtIteration.getModelMatrix();
        SI = Y(1, i) * Y(3, i);  % S * Ia
        D = Y(6, i);
        
        beta = paramsFit.beta;
        F = zeros(7, 1);  % vector just for SI
        F(1:2, 1) = [-beta / (Npop - D); beta / (Npop - D)].*SI;
        Y(:, i + 1) = RK4(iteration, Y(:, i), A, F, dt);
    end
end