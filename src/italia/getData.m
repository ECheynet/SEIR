function tableCOVIDItaly = getData(fileURL, dataLines)
% Collects the updated data of the COVID-19 pandemic in Italy from the
% Italian governement (https://github.com/pcm-dpc/COVID-19)

%% Input handling
% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
    if nargin < 1
        fileURL = 'https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni.csv';
    end
end

%% Setup the Import Options and import the data
% Set the number of columns
opts = delimitedTextImportOptions("NumVariables", 16);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types

%  opts.VariableNames = ["HospitalizedWithSymptoms", "HospitalizedInIntensiveCare", "Hospitalized",  	"HomeConfinement", 	"confirmed" ,	"activeCases", 	"newCases", "Recovered", 	"Deaths" ,	"totalCases"]
opts.VariableNames = ["Date"    , "CountryCode", "RegionCode", "RegionName", "Lat"   , "Long"  , "HospitalizedWithSymptoms", "HospitalizedInIntensiveCare", "Hospitalized", "HomeConfinement", "Quarantined", "variationQuarantined", "NewQuarantined","Recovered", "Deaths", "Confirmed" , "Swabs" ];
opts.VariableTypes = ["string", "string"     , "uint8"     , "string"    , "double", "double", "double"                  , "double"                     , "double"      , "double"         , "double"     , "double"        , "double"        , "double"   , "double", "double"    , "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
% Can create an error with older Matlab versions
% opts = setvaropts(opts, "Date", "InputFormat", "yyyy-MM-dd HH:mm:ss"); % Use "Format","yyyy-MM-dd'T'HH:mm:ss" from Mar 23, 2020 on

% Download the CSV file
websave('dummy.csv',fileURL);

% Import the data
fid = fopen('dummy.csv');
tableCOVIDItaly = readtable('dummy.csv', opts);
fclose(fid);
delete('dummy.csv')
end
