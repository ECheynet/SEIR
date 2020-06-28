function [tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID_US()
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

%% Number of days of data
Ndays = floor(datenum(now))-datenum(2020,01,22)-1; % minus one day because the data are updated with a delay of 24 h
address = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/';
ext = '.csv';
%% Options and names for confirmed
opts = delimitedTextImportOptions("NumVariables", Ndays+11);
opts.VariableNames = ["UID",           "iso2" ,       "iso3" ,   "code3" ,   "FIPS" ,   "Admin2" ,   "Province_State",    "Country_Region" ,   "Lat" ,    "Long_" ,   "Combined_Key", repmat("day",1,Ndays+1)];
opts.VariableTypes = [repmat("string",1,11), repmat("double",1,Ndays+1)];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";



filename = ['time_series_covid19_confirmed_US'];
fullName = [address,filename,ext];
urlwrite(fullName,'dummy.csv');
tableConfirmed =readtable('dummy.csv', opts);
delete('dummy.csv')
%% Options and names for deceased
%  One more row is used for the population!
%  Inconsistent format used by John Hopkins university

clear opts
opts = delimitedTextImportOptions("NumVariables", Ndays+12);
opts.VariableNames = ["UID",           "iso2" ,       "iso3" ,   "code3" ,   "FIPS" ,   "Admin2" ,   "Province_State",    "Country_Region" ,   "Lat" ,    "Long_" ,   "Combined_Key", "Population", repmat("day",1,Ndays+1)];
opts.VariableTypes = [repmat("string",1,11), repmat("double",1,Ndays+2)];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";


filename = ['time_series_covid19_deaths_US'];
fullName = [address,filename,ext];
urlwrite(fullName,'dummy.csv');
tableDeaths =readtable('dummy.csv', opts);
delete('dummy.csv')

%% Get time
time = datetime(2020,01,22):days(1):datetime(datestr(floor(datenum(now))), 'Locale', 'en_US')-datenum(1);
%% So far no data on recovered

tableRecovered = [];

end