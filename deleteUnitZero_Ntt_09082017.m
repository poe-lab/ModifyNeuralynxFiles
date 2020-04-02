function deleteUnitZero_Ntt_09082017
%% Select the Neuralynx tetrode file (.ntt):
working_dir=pwd;
dataFolder = [];
fileName = [];
fileSelectedCheck = 0;
while isequal(fileSelectedCheck,0)
    [fileName, dataFolder] = uigetfile({'*.ntt'}, 'Select the Neuralynx tetrode file');
    if isempty(fileName) || isempty(dataFolder)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        fileSelectedCheck = 1;
    end 
end
neuralynxFile = fullfile(dataFolder,fileName); %Full file path for .ntt file to be loaded
cd(working_dir);

%% Load the .ntt file:    
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
    Nlx2MatSpike(neuralynxFile, [1 1 1 1 1], 1, 1, [] );

%% Remove Unit 0 (unsorted) from the data:
cellNum = 0;
targetIdx = CellNumbers ~= cellNum;
Timestamps = Timestamps(targetIdx);
ScNumbers = ScNumbers(targetIdx);
CellNumbers = CellNumbers(targetIdx);
Features = Features(:, targetIdx);
Samples = Samples(:, :, targetIdx);

%% Write data to a new .ntt file:
newNttFile = fullfile(dataFolder, ['deleteUnitZero_' fileName]); %Full file path for new .ntt file w/o Unit 0
Mat2NlxSpike(newNttFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers,...
    CellNumbers, Features, Samples, Header);
end   
