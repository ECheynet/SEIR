function [A] = getA(gamma,delta,lambda,kappa)
    A = zeros(6);
    A(1,1) = 0;  % S
    A(2,2) = -gamma;  % E
    A(3,2:3) = [gamma,-delta];  % I
    A(4,3:4) = [delta,-kappa-lambda];  % Q
    A(5,4) = lambda;  % R
    A(6,4) = kappa;  % D
end