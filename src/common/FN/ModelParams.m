classdef ModelParams
   properties
      alpha
      beta
      gamma
      delta
      lambda
      kappa
      tau
      rho
      n
   end
   methods        
        function obj = ModelParams(coefficients)
            obj.alpha = coefficients(1);
            obj.beta = coefficients(2);
            obj.gamma = coefficients(3);
            obj.delta = coefficients(4);
            obj.lambda = coefficients(5);
            obj.kappa = coefficients(6);
            obj.tau = coefficients(7);
            obj.rho = coefficients(8);
            
            obj.n = length(coefficients);
        end
        function r = getAsVector(obj)
             r = [obj.alpha, ...
                obj.beta, ...
                obj.gamma, ...
                obj.delta, ...
                obj.lambda, ...
                obj.kappa, ...
                obj.tau, ...
                obj.rho ...
             ];
        end
        function lowerBounds = getLowerBounds(obj)
            lowerBounds = zeros(1, obj.n);
            
            lowerBounds(2) = 0.4;  % susceptible -> exposed
            lowerBounds(3) = 1/20; % 1 / latent period
        end
        function upperBounds = getUpperBounds(obj)
            upperBounds = ones(1, obj.n);
            
            upperBounds(3) = 1/15; % 1 / latent period
        end
        function A = getModelMatrix(obj)
            A = 1 % todo
        end
    end
end