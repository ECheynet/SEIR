function [A] = getModelMatrix(interactions)
    alpha = interactions(1); % unpack interactions
    gamma = interactions(2);
    delta = interactions(3);
    lambda = interactions(4);
    kappa = interactions(5);
    tau = interactions(6);
    rho = interactions(7);

    A = zeros(7);
    A(1, 1) = -alpha;  % out of  S
    A(2, 2) = -gamma;  % out of E
    A(3, 2:3) = [gamma, -delta-tau-rho];  % asym
    A(4, 3:4) = [delta, -kappa-lambda];  % test positive
    A(5, 4) = lambda;  % recover: lambda * Iq
    A(6, 3:4) = [rho, kappa];  % dead
    A(7, 1:3) = [alpha, 0, tau];  % not susciptable
end