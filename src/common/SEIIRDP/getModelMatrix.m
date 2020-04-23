function [A] = getModelMatrix(alpha, gamma, delta, lambda, kappa, tau, rho)
    A = zeros(7);
    A(1, 1) = -alpha;  % S
    A(2, 2) = -gamma;  % E
    A(3, 2:4) = [gamma, -lambda - kappa, -delta];  % Ia
    A(4, 3:4) = [delta, -tau - rho];  % Iq
    A(5, 3:4) = [lambda, tau];  % R
    A(6, 3:4) = [kappa, rho];  % D
    A(7, 1) = alpha;  % D
end