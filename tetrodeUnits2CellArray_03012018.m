function tetrodeUnits2CellArray_03012018
% This function imports all of the tetrode files in a folder and creates a
% cell array of units.
%--INPUT: NTT files in the selected folder
%--OUTPUT: .MAT file containing: Tetrode #, Unit # of tetrode, time of each
%spike

%% Load NTT files
% Select folder and get list of NTT files:
fileType = '*.ntt';
[dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
clear fileType

%% For each tetrode, load in the data and combine all spikes into one file:
allSpikeTS = [];
allID = [];
for i = 1:numberOfDataFiles
    fileName = strtrim(fileList(i,:)); %Removes any white space at end of file name string.
    tetrodeFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    
    %% Import the data:
    % Set up variable for sorted Neuralynx NTT file:
    [spikeTimes, tetrodeNum, cellNumber] = Nlx2MatSpike(tetrodeFile, [1 1 1 0 0], 0, 1, []); % Load only time stamps, cell #, and amplitudes
    
    %% Remove unsorted spikes:
    nonZerosIndex = find(cellNumber);       % Identify spikes with a unit assignment
    cellNumber = cellNumber(nonZerosIndex)'; % Remove cell # unsorted spikes
    spikeTimes = spikeTimes(nonZerosIndex)' ./ 1000000;   % Remove time stamps of unsorted spikes and convert to seconds
    tetrodeNum = tetrodeNum(nonZerosIndex)' + 1; %TT # starts at 0 in the file data
    cellIdentity = [tetrodeNum cellNumber];
    clear nonZerosIndex cellNumber tetrodeNum
    allSpikeTS = [allSpikeTS; spikeTimes];
    allID = [allID; cellIdentity];
    clear spikeTimes cellIdentity fileName
end
clear numberOFDataFiles

%% Sort all spike data by time stamp
[allSpikeTS, IX] = sort(allSpikeTS,1);
allID = allID(IX,:);
clear IX

%% Convert spike data to cell arrays
unit_ID = unique(allID, 'rows');
numUnits = size(unit_ID,1);
cellsOfUnits = cell(numUnits,1);
for i = 1:numUnits
    targetIdx = ismember(allID, unit_ID(i,:), 'rows');
    cellsOfUnits{i,1} = allSpikeTS(targetIdx);
    clear targetIdx
end
clear numUnits allID allSpikeTS

%% SAVE DATA
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
resultsFolder = dataFolder; %'Z:\Data Analysis\Optogenetic_AnalzedData\VLMC_Analyses';
save(fullfile(resultsFolder,['tetrodeUnitsTS', filename, '.mat']),...
    'fileList', 'unit_ID', 'cellsOfUnits');
