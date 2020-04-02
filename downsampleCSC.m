function downsampleCSC
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );
[m1,n1]=size(Samples);
newM = m1*n1;
Samples = reshape(Samples, newM, 1);
Samples = Samples(1:8:newM);
[m2,n2]=size(Samples);
newEnd = floor(m2/512);
shortEnd = newEnd * 512;
Samples = Samples(1:shortEnd);
Samples = reshape(Samples, m1, newEnd);
Timestamps = Timestamps(1:8:n1);
ChannelNumbers =ChannelNumbers(1:8:n1);
SampleFrequencies = SampleFrequencies(1:8:n1);
NumberOfValidSamples =NumberOfValidSamples(1:8:n1);
Timestamps = Timestamps(1:newEnd);
ChannelNumbers =ChannelNumbers(1:newEnd);
SampleFrequencies = SampleFrequencies(1:newEnd)./8;
NumberOfValidSamples =NumberOfValidSamples(1:newEnd);
% UNCOMMENT for NEW Neuralynx system files:
Header{13,1} = ['-SamplingFrequency ' num2str(SampleFrequencies(1))];
% UNCOMMENT for OLD Neuralynx system files:
% Header{14,1} = ['-SamplingFrequency ' num2str(SampleFrequencies(1))];

fixedNttFilename = ['downsampledTo' num2str(SampleFrequencies(1)) 'Hz_' CSCFilename ];
cscFile = fullfile(CSCFilePath, fixedNttFilename);
Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
    SampleFrequencies, NumberOfValidSamples, Samples, Header);
clear all
