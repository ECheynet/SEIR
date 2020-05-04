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
opts.VariableNames = ["Date", "_", "_", "_", "_", "_", "hospitalizedSymptomatic", "ICU", "hospitalized", "quarantined", "totPositive", "_", "_", "recovered", "dead", "totCases", "test"];
opts.VariableTypes = ["string", "string" , "uint8" , "string" , "double", "double", "double" , "double" , "double" , "double" , "double" , "double" , "double" , "double" , "double", "double" , "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.PreserveVariableNames = true;

% Specify variable properties
% Can create an error with older Matlab versions
% opts = setvaropts(opts, "Date", "InputFormat", "yyyy-MM-dd HH:mm:ss"); % Use "Format","yyyy-MM-dd'T'HH:mm:ss" from Mar 23, 2020 on

% Download the CSV file
% uncomment to to update data
websave('data.csv', fileURL);

% Import the data
fid = fopen('data.csv');
tableCOVIDItaly = readtable('data.csv', opts);
fclose(fid);
% delete('dummy.csv')
end
