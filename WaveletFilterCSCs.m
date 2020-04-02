function WaveletFilterCSCs
% Select folder where Neuralynx files are located:
fileSelectedCheck = 0;
while isequal(fileSelectedCheck,0)
    fileType = '*.ncs*';
    [dataFolder, fileList, numberOfDataFiles] = batchLoadFiles(fileType);
    if isempty(dataFolder) || isempty(fileList)
        uiwait(errordlg('You need to select a folder with data. Please try again',...
            'ERROR','modal'));
    else
        fileSelectedCheck = 1;
    end 
end    
mkdir(dataFolder, 'Wavelet_Filtered')
% Define the high pass cut off frequency of the wavelet.
% The cutoff frequency = samplingrate/(2^(maxlevel+1)):
maxlevel = 6;

for m = 1:numberOfDataFiles
    fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    waveletFile = fullfile(dataFolder, 'Wavelet_Filtered', ['waveFilt_' fileName]); %Full file path for new data file with reset time stamps
    
    % Load the file to apply the wavelet filter to the signal:
    [Timestamps, ChannelNumbers, SampleFreq, numSamples, Samples, Header] =...
        Nlx2MatCSC(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
    
    % Reshape data into a column vector:
    [m1,n1] = size(Samples);
    newM = m1*n1;
    Samples = reshape(Samples, newM, 1);
    
    % Apply wavelet filter:
    Samples = wavefilter(Samples', maxlevel);
    
    %Convert back to NLX format:
    Samples = reshape(Samples, m1, n1); 
    
    % Save wavelet filtered data as a new Neuralynx CSC (.NCS) file:
    Mat2NlxCSC(waveletFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
        SampleFreq, numSamples, Samples, Header);
    clear Timestamps ChannelNumbers SampleFreq numSamples Samples Header
end