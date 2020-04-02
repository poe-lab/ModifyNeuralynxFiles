function CommonAverageRef_DiffShankCSCs07032017
% Select folder where Neuralynx files are located to calculate CAR:
fileSelectedCheck = 0;
dialogTitle = 'Select folder with CSCs to be used to calculate CAR';
while isequal(fileSelectedCheck,0)
    fileType = '*.ncs*';
    [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType,dialogTitle);
    if isempty(dataFolder) || isempty(fileList)
        uiwait(errordlg('You need to select a folder with data. Please try again',...
            'ERROR','modal'));
    else
        fileSelectedCheck = 1;
    end 
end    

sumCscs = [];
for m = 1:numberOfDataFiles
    fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    % load the samples to add to sum:    
    Samples = Nlx2MatCSC(neuralynxFile, [0 0 0 0 1], 0, 1, []);
    if isequal(m, 1)
        sumCscs = Samples;
    else
        sumCscs = sumCscs + Samples;
    end
    clear Samples fileName neuralynxFile
end

% Calculate the common average reference signal:
CAR = sumCscs./numberOfDataFiles;
carDataFolder = dataFolder;
clear sumCscs numberOfDataFiles dataFolder fileList

% Select folder where Neuralynx files are located to re-reference to CAR:
fileSelectedCheck = 0;
dialogTitle = 'Select folder with CSCs to be referenced to CAR';
while isequal(fileSelectedCheck,0)
    fileType = '*.ncs*';
    [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType,dialogTitle);
    if isempty(dataFolder) || isempty(fileList)
        uiwait(errordlg('You need to select a folder with data. Please try again',...
            'ERROR','modal'));
    elseif isequal(carDataFolder, dataFolder)
        uiwait(errordlg('You need to select a different for than the CAR data folder. Please try again',...
            'ERROR','modal'));
    else
        fileSelectedCheck = 1;
    end 
end

mkdir(dataFolder, 'CAR') % creates new subdirectory for CAR referenced CSCs

for m = 1:numberOfDataFiles
    fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    carFile = fullfile(dataFolder, 'CAR', ['CAR_' fileName]); %Full file path for new data file with reset time stamps
    % load the file and create a new version referenced to common average:
    [Timestamps, ChannelNumbers, SampleFreq, numSamples, Samples, Header] =...
        Nlx2MatCSC(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
    Samples = Samples - CAR;
    Mat2NlxCSC(carFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
        SampleFreq, numSamples, Samples, Header);
    clear Timestamps ChannelNumbers SampleFreq numSamples Samples Header
end