[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );
fixedNcsFilename = ['Part1' CSCFilename ];
cscFile = fullfile(CSCFilePath, fixedNcsFilename);
Pt1NumberOfValidSamples = NumberOfValidSamples(1:180992);
Pt2NumberOfValidSamples = NumberOfValidSamples(180993:end);
Pt1SampleFrequencies = SampleFrequencies(1:180992);
Pt2SampleFrequencies = SampleFrequencies(180993:end);
Pt1Samples = Samples(:,1:180992);
Pt2Samples = Samples(:,180993:end);
Pt1Timestamps = Timestamps(1:180992);
Pt2Timestamps = Timestamps(180993:end);
Pt1ChannelNumbers = ChannelNumbers(1:180992);
Pt2ChannelNumbers = ChannelNumbers(180993:end);
Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Pt1Timestamps, Pt1ChannelNumbers,...
Pt1SampleFrequencies, Pt1NumberOfValidSamples, Pt1Samples, Header);
fixedNcsFilename = ['Part2' CSCFilename ];
cscFile = fullfile(CSCFilePath, fixedNcsFilename);
Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Pt2Timestamps, Pt2ChannelNumbers,...
Pt2SampleFrequencies, Pt2NumberOfValidSamples, Pt2Samples, Header);