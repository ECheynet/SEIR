function [A] = getModelMatrix(interactions)
    alpha = interactions(1); % unpack interactions
    gamma = interactions(2);
    delta = interactions(3);
    lambda = interactions(4);
    kappa = interactions(5);
    tau = interactions(6);

    A = zeros(7);
    A(1, 1) = -alpha;  % out of  S
    A(2, 2) = -gamma;  % out of E
    A(3, 2:3) = [gamma, -delta-tau];  % asym: + gamma * E - delta
    A(4, 3:4) = [delta, -kappa-lambda];  % test positive: + delta * asym - lambda
    A(5, 3:4) = [tau, lambda];  % recover: +tau * asym + lambda * Iq
    A(6, 4) = kappa;  % dead: + kappa * Iq
    A(7, 1) = alpha;  % not susciptable: + alpha * S
end