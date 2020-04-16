function [rmse, nrmse] = mof(realData, yourPrediction)
% Measures of error

%% parse
n = min(length(realData), length(yourPrediction));
realData = realData(1:n);
yourPrediction = yourPrediction(1:n);

%% RMSE
rmse = sqrt(sum(realData(:) - yourPrediction(:))^2/n);

%% Normalized RMSE
nrmse = rmse / mean(realData);
end

