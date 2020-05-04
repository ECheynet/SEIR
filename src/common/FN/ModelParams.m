classdef ModelParams
   properties
      alpha % protection rate
      beta % S -> E (by coming in contact with asymp)
      gamma % 1 / latent time
      delta % asymp -> test positive
      lambda % recovery rate (when being symptomatic)
      kappa % death rate (when being symptomatic)
      tau % asym -> recover
      rho % death rate (when being asymptomatic)
      n  % number of coefficients
   end
   methods        
        function obj = ModelParams(coefficients)
            obj.alpha = coefficients(1);
            obj.beta = coefficients(2);
            obj.gamma = coefficients(3);
            obj.delta = coefficients(4);
            obj.lambda = coefficients(5:6);
            obj.kappa = coefficients(7);
            obj.tau = coefficients(8);
            obj.rho = coefficients(9);
            
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
            A = zeros(7);
            
            A(1, 1) = -obj.alpha;  % out of  S
            A(2, 2) = -obj.gamma;  % out of E
            A(3, 2:3) = [obj.gamma, -obj.delta-obj.tau-obj.rho];  % asym
            A(4, 3:4) = [obj.delta, -obj.kappa-obj.lambda];  % test positive
            A(5, 4) = obj.lambda;  % recover: lambda * Iq
            A(6, 3:4) = [obj.rho, obj.kappa];  % dead
            A(7, 1:3) = [obj.alpha, 0, obj.tau];  % not susciptable
        end
    end
end