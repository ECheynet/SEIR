function [A] = getModelMatrix(alpha, gamma, delta, lambda, kappa)
    A = zeros(7);
    A(1, 1) = -alpha;  % out of  S
    A(2, 2) = -gamma;  % out of E
    A(3, 2:3) = [gamma, -delta];  % + gamma * E - delta
    A(4, 3:4) = [delta, -kappa-lambda];  % + delta * Ia - lambda
    A(5, 4) = lambda;  % + lambda * Iq
    A(6, 4) = kappa;  % + kappa * Iq
    A(7, 1) = alpha;  % + alpha * S
end