function JoinCscFiles
working_dir = pwd;
% Load first CSC file:

[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick  1st CSC file.'},'Select 1st CSC Data File');
cd(working_dir);
CSCFile = fullfile(CSCFilePath, CSCFilename);

[Timestamps1, ChannelNumbers1, SampleFrequencies1,NumberOfValidSamples1, Samples1, Header1]...
    = Nlx2MatCSC(CSCFile,[1 1 1 1 1], 1, 1, [] );

% Load 2nd CSC file's header:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick  2nd CSC file.'},'Select 2nd CSC Data File');
CSCFile = fullfile(CSCFilePath, CSCFilename);

[Timestamps2, ChannelNumbers2, SampleFrequencies2,NumberOfValidSamples2, Samples2, Header2]...
    = Nlx2MatCSC(CSCFile,[1 1 1 1 1], 1, 1, [] );

% Create header for joined CSC file:
Header1{4,1} =  Header2{4,1};

fixedCSCFilename = ['Joined' CSCFilename];
CSCFile = fullfile(CSCFilePath, fixedCSCFilename);
Mat2NlxCSC(CSCFile, 0, 1, 1, [1 1 1 1 1 1], [Timestamps1 Timestamps2], [ChannelNumbers1 ChannelNumbers2],...
    [SampleFrequencies1 SampleFrequencies2], [NumberOfValidSamples1 NumberOfValidSamples2], [Samples1 Samples2], Header1);



