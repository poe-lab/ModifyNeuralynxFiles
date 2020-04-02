function filterCSC
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );
[m1,n1]=size(Samples);
newM = m1*n1;
Samples = reshape(Samples, newM, 1);
%Now filter signal:
[z, p, k] = ellip(7,1,60, [1 30]/(SampleFrequencies(1)/2),'bandpass');
[sos, g] = zp2sos(z,p,k);
Samples = filtfilt(sos,g, Samples);

%Convert back to NLX format:
Samples = reshape(Samples, m1, n1);



% UNCOMMENT for NEW Neuralynx system files:
%Header{13,1} = ['-SamplingFrequency ' num2str(SampleFrequencies(1))];
% UNCOMMENT for OLD Neuralynx system files:
% Header{14,1} = ['-SamplingFrequency ' num2str(SampleFrequencies(1))];

filteredCSCFilename = ['1to30Hz_Filtered' CSCFilename ];
cscFile = fullfile(CSCFilePath, filteredCSCFilename);
Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
    SampleFrequencies, NumberOfValidSamples, Samples, Header);
clear all
