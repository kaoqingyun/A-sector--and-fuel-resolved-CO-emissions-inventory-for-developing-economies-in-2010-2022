%% Six-factor LMDI decomposition, 2010-2022
% The change in country-level CO2 emissions is decomposed into:
% population, GDP per capita, industrial structure, energy intensity,
% energy mix, and emission factor effects.
%
% Required variables after column-name standardization:
% Country, Year, Sector, E, C, G_sec, G_total, and P.
%
% The observations for the start and end years are assumed to have the
% same sector-energy coverage and row order.

clc;
clear;
close all;

%% Analysis period
startYear = 2010;
endYear   = 2022;

%% Import data
fprintf('Reading input data...\n');
[fileName, filePath] = uigetfile( ...
    {'*.csv;*.xlsx', 'Data files (*.csv, *.xlsx)'}, ...
    'Select the input data file');

if isequal(fileName, 0)
    fprintf('File selection cancelled.\n');
    return;
end

T = readtable(fullfile(filePath, fileName));

%% Standardize variable names
nameMap = {
    'Gi',      'G_sec';
    '产值',    'G_sec';
    'G_',      'G_total';
    '总GDP',   'G_total';
    'G (',     'G_total';
    'Eij',     'E';
    '能耗',    'E';
    'Cij',     'C';
    '排放',    'C';
    'Pop',     'P';
    '人口',    'P';
    'P_',      'P';
    'P (',     'P';
    'Year',    'Year';
    'Country', 'Country';
    'Sector',  'Sector'
};

originalNames = T.Properties.VariableNames;
renamed = false(size(originalNames));
standardNames = {
    'Country', 'Year', 'Sector', 'E', 'C', 'G_sec', 'G_total', 'P'
};

for i = 1:numel(originalNames)
    renamed(i) = any(strcmp(originalNames{i}, standardNames));
end

for j = 1:size(nameMap, 1)
    sourceKey = nameMap{j, 1};
    targetName = nameMap{j, 2};

    for i = 1:numel(originalNames)
        if ~renamed(i) && contains(originalNames{i}, sourceKey, ...
                'IgnoreCase', true)
            T.Properties.VariableNames{i} = targetName;
            renamed(i) = true;
        end
    end
end

if ~ismember('P', T.Properties.VariableNames)
    for i = 1:numel(originalNames)
        if ~renamed(i) && startsWith(originalNames{i}, 'P', ...
                'IgnoreCase', true)
            T.Properties.VariableNames{i} = 'P';
            break;
        end
    end
end

requiredVariables = {
    'Country', 'Year', 'Sector', 'E', 'C', 'G_sec', 'G_total', 'P'
};
missingVariables = setdiff(requiredVariables, T.Properties.VariableNames);

if ~isempty(missingVariables)
    error('Missing required variables: %s', ...
        strjoin(missingVariables, ', '));
end

%% Convert numeric variables
numericVariables = {'E', 'C', 'G_sec', 'G_total', 'P', 'Year'};

for i = 1:numel(numericVariables)
    variableName = numericVariables{i};
    rawData = T.(variableName);

    if iscell(rawData) || isstring(rawData)
        cleanData = zeros(height(T), 1);

        for r = 1:height(T)
            valueText = char(string(rawData(r)));
            valueText = strrep(valueText, ',', '');
            valueText = strrep(valueText, ' ', '');

            if isempty(valueText) || strcmpi(valueText, '-') || ...
                    strcmpi(valueText, 'null')
                valueText = '0';
            end

            value = str2double(valueText);
            if isnan(value)
                value = 0;
            end
            cleanData(r) = value;
        end

        T.(variableName) = cleanData;
    else
        T.(variableName)(isnan(T.(variableName))) = 0;
    end
end

% A small positive value is used to permit logarithmic calculations.
epsilon = 1e-10;
T.E(T.E == 0) = epsilon;
T.C(T.C == 0) = epsilon;

T.Country = string(T.Country);
T.Sector = string(T.Sector);

%% Country-level LMDI decomposition
fprintf('Calculating the decomposition for %d-%d...\n', ...
    startYear, endYear);

countryList = unique(T.Country, 'stable');
resultData = cell(numel(countryList), 11);
excludedCountries = strings(numel(countryList), 1);

resultCount = 0;
excludedCount = 0;

logMean = @(x, y) (x - y) ./ (log(x) - log(y));

for i = 1:numel(countryList)
    countryName = countryList(i);
    countryData = T(T.Country == countryName, :);

    dataStart = countryData(countryData.Year == startYear, :);
    dataEnd = countryData(countryData.Year == endYear, :);

    if isempty(dataStart) || isempty(dataEnd)
        excludedCount = excludedCount + 1;
        excludedCountries(excludedCount) = countryName;
        continue;
    end

    invalidKeyData = ...
        any(dataStart.G_total == 0) || any(dataEnd.G_total == 0) || ...
        any(dataStart.P == 0)       || any(dataEnd.P == 0)       || ...
        any(dataStart.G_sec == 0)   || any(dataEnd.G_sec == 0);

    if invalidKeyData
        excludedCount = excludedCount + 1;
        excludedCountries(excludedCount) = countryName;
        continue;
    end

    populationEnd = dataEnd.P;
    populationStart = dataStart.P;
    gdpEnd = dataEnd.G_total;
    gdpStart = dataStart.G_total;
    sectorOutputEnd = dataEnd.G_sec;
    sectorOutputStart = dataStart.G_sec;
    energyEnd = dataEnd.E;
    energyStart = dataStart.E;
    emissionsEnd = dataEnd.C;
    emissionsStart = dataStart.C;

    weights = logMean(emissionsEnd, emissionsStart);
    equalEmissions = emissionsEnd == emissionsStart;
    weights(equalEmissions) = emissionsEnd(equalEmissions);

    populationEffect = sum(weights .* ...
        log(populationEnd ./ populationStart));

    activityEffect = sum(weights .* log( ...
        (gdpEnd ./ populationEnd) ./ ...
        (gdpStart ./ populationStart)));

    structureEffect = sum(weights .* log( ...
        (sectorOutputEnd ./ gdpEnd) ./ ...
        (sectorOutputStart ./ gdpStart)));

    sectorList = unique(dataEnd.Sector, 'stable');
    sectorEnergyEnd = zeros(size(energyEnd));
    sectorEnergyStart = zeros(size(energyStart));

    for k = 1:numel(sectorList)
        sectorMask = dataEnd.Sector == sectorList(k);
        sectorEnergyEnd(sectorMask) = sum(energyEnd(sectorMask));
        sectorEnergyStart(sectorMask) = sum(energyStart(sectorMask));
    end

    intensityEffect = sum(weights .* log( ...
        (sectorEnergyEnd ./ sectorOutputEnd) ./ ...
        (sectorEnergyStart ./ sectorOutputStart)));

    energyMixEffect = sum(weights .* log( ...
        (energyEnd ./ sectorEnergyEnd) ./ ...
        (energyStart ./ sectorEnergyStart)));

    emissionFactorEffect = sum(weights .* log( ...
        (emissionsEnd ./ energyEnd) ./ ...
        (emissionsStart ./ energyStart)));

    calculatedChange = populationEffect + activityEffect + ...
        structureEffect + intensityEffect + energyMixEffect + ...
        emissionFactorEffect;

    observedChange = sum(emissionsEnd - emissionsStart);
    residual = calculatedChange - observedChange;

    resultCount = resultCount + 1;
    resultData(resultCount, :) = {
        char(countryName), endYear, populationEffect, activityEffect, ...
        structureEffect, intensityEffect, energyMixEffect, ...
        emissionFactorEffect, calculatedChange, observedChange, residual
    };
end

resultData = resultData(1:resultCount, :);
excludedCountries = excludedCountries(1:excludedCount);

%% Report and export results
fprintf('\nDecomposition summary\n');
fprintf('Analysis period: %d-%d\n', startYear, endYear);
fprintf('Countries included: %d\n', resultCount);
fprintf('Countries excluded: %d\n', excludedCount);

if excludedCount > 0
    fprintf('Excluded countries:\n');
    disp(excludedCountries);
end

if resultCount > 0
    resultTable = cell2table(resultData, 'VariableNames', {
        'Country', 'Year', 'Pop', 'GDP_Per_Capita', 'Structure', ...
        'Intensity', 'Mix', 'EmissionFactor', 'Calc_Sum', ...
        'Real_Sum', 'Residual'
    });

    disp(resultTable(1:min(5, height(resultTable)), :));

    outputFile = sprintf('LMDI_%d_%d_Results.csv', ...
        startYear, endYear);
    writetable(resultTable, outputFile);

    fprintf('Results saved to %s\n', outputFile);
else
    fprintf('No country met the calculation requirements.\n');
end
