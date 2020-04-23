function [A] = getIterationMatrix(gamma, delta, lambda, kappa, tau, rho)
    A = zeros(6);
    A(1, 1) = 0;  % S
    A(2, 2) = [-gamma - delta];  % E
    A(3, 2:3) = [gamma, -lambda - kappa];  % Ia
    A(4, 2:4) = [delta, 0, -tau - rho];  % Iq
    A(5, 3:4) = [lambda, tau];  % R
    A(6, 3:4) = [kappa, rho];  % D
end