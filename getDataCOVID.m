function [tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID()
% The function [tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID
% collect the updated data from the COVID-19 epidemy from the
% John Hopkins university [1]
% 
% References:
% [1] https://github.com/CSSEGISandData/COVID-19
% 
% Author: E. Cheynet - Last modified - 20-03-2020
% 
% see also fit_SEIQRDP.m SEIQRDP.m
%% Options and names
Ndays = floor(datenum(now))-datenum(2020,01,22)-1; % minus one day because the data are updated with a delay of 24 h
opts = delimitedTextImportOptions("NumVariables", Ndays+5);
opts.VariableNames = ["ProvinceState", "CountryRegion", "Lat", "Long", repmat("data",1,Ndays+1)];
opts.VariableTypes = ["string", "string", repmat("double",1,Ndays+3)];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
% opts = setvaropts(opts, ["ProvinceState", "CountryRegion"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["ProvinceState", "CountryRegion"], "EmptyFieldRule", "auto");
%% Import the data
status = {'confirmed','deaths','recovered'};
address = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/';
ext = '.csv';
for ii=1:numel(status)
    
    filename = ['time_series_covid19_',status{ii},'_global'];
    fullName = [address,filename,ext];
%     disp(fullName)
    urlwrite(fullName,'dummy.csv');
    
    if strcmpi(status{ii},'Confirmed')
        tableConfirmed =readtable('dummy.csv', opts);
        
    elseif strcmpi(status{ii},'Deaths')
        tableDeaths =readtable('dummy.csv', opts);
        
    elseif strcmpi(status{ii},'Recovered')
        tableRecovered =readtable('dummy.csv', opts);
    else
        error('Unknown status')
    end
end

time = datetime(2020,01,22):days(1):datetime(datestr(floor(datenum(now))))-datenum(1);

% fid = fopen('dummy.csv');
% time = textscan(fid,repmat('%s',1,size(tableConfirmed,2)), 1, 'Delimiter',',');
% time(1:12)=[];
% time = datetime([time{1:end}])+years(2000);
% fclose(fid);
% delete('dummy.csv')

end