function fixHeaderCSC(batchProcess)
% Good CSC file header:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
goodHeader = Nlx2MatCSC(cscFile, [0 0 0 0 0], 1, 1, [] );

%% Load CSC files that need the header fixed:
working_dir=pwd;
if batchProcess
    % Select folder and get list of CSC files with bad header (i.e., system crashed during recording):
    fileType = '*.ncs';
    [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
else
    dataFolder = [];
    fileName = [];
    fileSelectedCheck = 0;
    % CSC file with bad header (i.e., system crashed during recording):
    while isequal(fileSelectedCheck,0)
        [CSCFilename, dataFolder] = uigetfile({'*.ncs',...
            'Pick CSC files.'},'Select Continuously Sampled Channel File');
        if isempty(fileName) || isempty(dataFolder)
            uiwait(errordlg('You need to select a file. Please try again',...
                'ERROR','modal'));
        else
            fileSelectedCheck = 1;
        end 
    end
    cd(working_dir);
    numberOfDataFiles = 1;
end

for i= 1:numberOfDataFiles
    if batchProcess
        CSCFilename = strtrim(fileList(i,:)); %Removes any whites space at end of file name string.
    end
    cscFile = fullfile(dataFolder,CSCFilename);
    [Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );

    % UNCOMMENT for NEW Neuralynx system files:
    % Header{13,1} = ['-SamplingFrequency ' num2str(SampleFrequencies(1))];
    % UNCOMMENT for OLD Neuralynx system files:
    Header{14,1} = goodHeader{14,1}; %Sampling Frequency
    Header{16,1} = goodHeader{16,1}; %AD Bit Volts
    Header{21,1} = goodHeader{21,1}; %Input range
    Header{24,1} = goodHeader{24,1}; %Amplitude Low Cut (system filter)
    Header{25,1} = goodHeader{25,1}; %Amplitude Hi Cut (system filter)
    Header{26,1} = goodHeader{26,1}; %Amplitude Gain
    fixedNcsFilename = ['fixedHdr' CSCFilename ];
    cscFile = fullfile(CSCFilePath, fixedNcsFilename);
    Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
        SampleFrequencies, NumberOfValidSamples, Samples, Header);
end
clear all
