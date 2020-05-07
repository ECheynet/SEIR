function [tableConfirmed,tableDeaths,tableRecovered,time] = getDataCOVID_FRA()
% The function [tableConfirmed,tableDeaths,tableRecovered,time] =
% getDataCOVID_FRA
% collect the updated data from the COVID-19 epidemy in France from [1]
% 
% References:
% [1]  https://github.com/cedricguadalupe/FRANCE-COVID-19 
% 
% Author: E. Cheynet - Last modified - 07-05-2020
% 
% see also fit_SEIQRDP.m SEIQRDP.m

%% Options and names

Nregions = 20;

opts = delimitedTextImportOptions('NumVariables', Nregions);
opts.VariableNames = ["Date", "AuvergneRhoneAlpes", "BourgogneFrancheComte", "Bretagne", "CentreValdeLoire", "Corse", "GrandEst", "HautsdeFrance", "IledeFrance", "Normandie", "NouvelleAquitaine", "Occitanie", "PaysdelaLoire", "ProvenceAlpesCotedAzur", "Guadeloupe", "Martinique", "Guyane", "LaRunion", "Mayotte",  "Total"];opts.VariableTypes(1) = {'string'};
opts.VariableTypes(1) = {'datetime'};
opts.VariableTypes(2:Nregions) = {'double'};
% Specify file level properties
opts.ExtraColumnsRule = 'ignore';
opts.EmptyLineRule = 'read';
% Specify variable properties
% opts = setvaropts(opts, ['ProvinceState', 'CountryRegion'], 'WhitespaceRule', 'preserve');
% opts = setvaropts(opts, ['ProvinceState', 'CountryRegion'], 'EmptyFieldRule', 'auto');


%% Import the data

status = {'confirmed','deaths','recovered'};
address = 'https://raw.githubusercontent.com/cedricguadalupe/FRANCE-COVID-19/master/';
ext = '.csv';


for ii=1:numel(status)
    
    filename = ['france_coronavirus_time_series-',status{ii}];
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


fid = fopen('dummy.csv');
time = tableRecovered.Date(1:end);
fclose(fid);

delete('dummy.csv')

end

