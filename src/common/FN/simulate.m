function [Y] = simulate(params, Y, Npop, t)
    N = numel(t);
    dt = median(diff(t));
    iteration = @(Y,A,F) A*Y + F;        
    A = getModelMatrix(params);  % ModelParams(params).
     
    for i=1:N-1
        SI = Y(1, i) * Y(3, i);  % S * Ia
        D = Y(6, i);
        
        F = zeros(7, 1);  % vector just for SI
        beta = params(2);  % params.beta
        F(1:2, 1) = [-beta / (Npop - D); beta / (Npop - D)].*SI;
        
        Y(:, i + 1) = RK4(iteration, Y(:, i), A, F, dt);
    end
end
