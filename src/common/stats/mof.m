function [rmse, nrmse] = mof(realData, yourPrediction)
% Measures of error

%% parse
n = min(length(realData), length(yourPrediction));
realData = realData(1:n);
yourPrediction = yourPrediction(1:n);

%% MSE
mse = mean((realData - yourPrediction).^2);

%% RMSE
rmse = sqrt(mse);

%% Normalized RMSE
nrmse = rmse / mean(realData);
end

