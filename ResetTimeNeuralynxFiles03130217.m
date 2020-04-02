function ResetTimeNeuralynxFiles03130217(batchProcess)
working_dir=pwd;
if batchProcess
    % Select folder where Neuralynx files are located:
    fileSelectedCheck = 0;
    while isequal(fileSelectedCheck,0)
        fileType = '*.n*';
        [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
        if isempty(dataFolder) || isempty(fileList)
            uiwait(errordlg('You need to select a folder with data. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end    
else
    numberOfDataFiles = 1;
    dataFolder = [];
    fileName = [];
    fileSelectedCheck = 0;
    
    % Select a single EDF file:
    while isequal(fileSelectedCheck,0)
        [fileName, dataFolder] = uigetfile({'*.ntt'; '*.ncs'; '*.nev'; '*.nvt'}, 'Select the Neuralynx data file');
        if isempty(fileName) || isempty(dataFolder)
            uiwait(errordlg('You need to select a file. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    cd(working_dir);
end

%% Find the start time stamp of each Neuralynx file:
startTime = NaN(numberOfDataFiles, 1);
for m = 1:numberOfDataFiles
    if batchProcess
        fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    end
    
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    if endsWith(neuralynxFile, '.ntt', 'IgnoreCase', true)
        startTime(m) = Nlx2MatSpike(neuralynxFile, [1 0 0 0 0], 0, 3, 1);
    elseif endsWith(neuralynxFile, '.ncs', 'IgnoreCase', true)
        startTime(m) = Nlx2MatCSC(neuralynxFile, [1 0 0 0 0], 0, 3, 1);
    elseif endsWith(neuralynxFile, '.nev', 'IgnoreCase', true)
        startTime(m) = Nlx2MatEV(neuralynxFile, [1 0 0 0 0], 0, 3, 1);
    elseif endsWith(neuralynxFile, '.nvt', 'IgnoreCase', true)
        startTime(m) = Nlx2MatVT(neuralynxFile, [1 0 0 0 0 0], 0, 3, 1);
    end    
end
% Find the earliest start time in the recording:
resetTime = min(startTime, [], 'omitnan') -1;
mkdir(dataFolder, 'ResetTime')

%% Create new files with the new start time:
for m = 1:numberOfDataFiles
    if batchProcess
        fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    end
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    resetFile = fullfile(dataFolder, 'ResetTime', ['tsReset' fileName]); %Full file path for new data file with reset time stamps
    % load the file and create a new version with reset time stamps:
    if endsWith(neuralynxFile, '.ntt', 'IgnoreCase', true)
        [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] =...
            Nlx2MatSpike(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
        Timestamps = Timestamps - resetTime;
        Mat2NlxSpike(resetFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers,...
            CellNumbers, Features, Samples, Header);
        clear Timestamps ScNumbers CellNumbers Features Samples Header
        
    elseif endsWith(neuralynxFile, '.ncs', 'IgnoreCase', true)
        [Timestamps, ChannelNumbers, SampleFreq, numSamples, Samples, Header] =...
            Nlx2MatCSC(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
        Timestamps = Timestamps - resetTime;
        Mat2NlxCSC(resetFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
            SampleFreq, numSamples, Samples, Header);
        clear Timestamps ChannelNumbers SampleFreq numSamples Samples Header
        
    elseif endsWith(neuralynxFile, '.nev', 'IgnoreCase', true)
        [Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] =...
            Nlx2MatEV(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
        Timestamps = Timestamps - resetTime;
        Mat2NlxEV(resetFile, 0, 1, [], [1 1 1 1 1],...
            Timestamps, EventIDs, TTLs, Extras, EventStrings, Header);
        clear Timestamps EventIDs TTLs Extras EventStrings Header
        
    elseif endsWith(neuralynxFile, '.nvt', 'IgnoreCase', true)
        [Timestamps, X, Y, Angles, Targets, Points, Header] =...
            Nlx2MatVT(neuralynxFile, [1 1 1 1 1 1], 1, 1, [] );
        Timestamps = Timestamps - resetTime;
        Mat2NlxVT(resetFile, 0, 1, [], [1 1 1 1 1 1],...
            Timestamps, X, Y, Angles, Targets, Points, Header);
        clear Timestamps X Y Angles Targets Points Header
        
    end    
end

end