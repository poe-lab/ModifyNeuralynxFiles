function downsampleClosestto1kHzCSC_v09172018
%% Load the Neuralynx CSC (.NCS) file:
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );
OriginalCSCFilename = CSCFilename;

%% Reduce variables to only necessary information
chanNum = ChannelNumbers(1);
clear ChannelNumbers
Fs = SampleFrequencies(1);
clear SampleFrequencies
numValidSamp = NumberOfValidSamples(1);
clear NumberOfValidSamples

%% Check for out of order time stamps due to memory buffering error:
modeDiffTS = [];
testOrderTS = 1;
while testOrderTS == 1    
    outOfOrderTS = find(diff(Timestamps) < 0, 1);
    if isempty(outOfOrderTS)
        testOrderTS = 0;
    else
        if isempty(modeDiffTS)
            modeDiffTS = mode(diff(Timestamps));
        end
        jumpToIdx = find(Timestamps >= (Timestamps(outOfOrderTS) + modeDiffTS), 1);
        Timestamps= [Timestamps(1:outOfOrderTS), Timestamps(jumpToIdx:end)];
        Samples = [Samples(:, 1:outOfOrderTS), Samples(:, jumpToIdx:end)];
    end
end

%% Check for breaks in the data:
A= min(diff(Timestamps));
T=find(diff(Timestamps)>1.1*A);
if ~isempty(T)
    [Samples, Timestamps] = fillTimeGaps(T, Samples, Timestamps);
    CSCFilename = ['FillGaps_'  CSCFilename];
end
clear T A

%% Convert the samples into a vector:
[m1,n1]=size(Samples);
newM = m1*n1;
Samples = reshape(Samples, newM, 1);

%% Find high pass filter in Header:
targ= strfind(Header,'-DspLowCutFrequency');
for i=1:length(targ)
    targIdx(i)= isempty(targ{i}); %#ok<AGROW>
end
LowFreqCutIdx = find(targIdx==0);   
clear targ targIdx
LowFreqCut = str2double(strrep(Header{LowFreqCutIdx,1}, '-DspLowCutFrequency ', ''));

%% Filter data if the current low cut filter < 0.5 Hz:
newLowFreqCut = 0.5;
if LowFreqCut < newLowFreqCut
    [z, p, k] = ellip(7,1,60, newLowFreqCut/(Fs/2),'high');
    [sos, g] = zp2sos(z,p,k);
    Samples = filtfilt(sos,g, Samples);
    % Update Low Cut frequency in Header:
    Header{LowFreqCutIdx,1} = ['-DspLowCutFrequency ' num2str(newLowFreqCut)];
    CSCFilename = ['LowCutPt5Hz_'  CSCFilename];
end

%% Calculate new sampling frequency closest to 1 kHz:
downSampFactor = floor(Fs/1000);
newFs = Fs/downSampFactor;

if downSampFactor == 1
    % If the sampling rate is as close to 1 kHz as possible, just convert
    % samples back to 512 x n for Neuralynx:
    [m2,~]=size(Samples);
    newEnd = floor(m2/512);
    shortEnd = newEnd * 512;
    Samples = Samples(1:shortEnd);
    Samples = reshape(Samples, m1, newEnd);
else        
    %% Find low pass filter in Header:
    targ= strfind(Header,'-DspHighCutFrequency');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    HighFreqCutIdx = find(targIdx==0);   
    clear targ targIdx
    HighFreqCut = str2double(strrep(Header{HighFreqCutIdx,1}, '-DspHighCutFrequency ', ''));
    
    %% Filter data if the current high cut filter *3 is > the down-sampled Fs:
    if newFs < HighFreqCut*3
        newHighFreqCut = 300;
        [z, p, k] = ellip(7,1,60, newHighFreqCut/(Fs/2),'low');
        [sos, g] = zp2sos(z,p,k);
        Samples = filtfilt(sos,g, Samples);
        % Update High Cut frequency in Header:
        Header{HighFreqCutIdx,1} = ['-DspHighCutFrequency ' num2str(newHighFreqCut)];
    end
    %% Down-sample the signal and time stamps closest to 1000 samples/second:
    Samples = Samples(1:downSampFactor:newM);
    [m2,~]=size(Samples);
    newEnd = floor(m2/512);
    shortEnd = newEnd * 512;
    Samples = Samples(1:shortEnd);
    Samples = reshape(Samples, m1, newEnd);
    Timestamps = Timestamps(1:downSampFactor:n1);

    Timestamps = Timestamps(1:newEnd);
    
    %% Update sampling frequency in Header:
    targ= strfind(Header,'-SamplingFreq');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    headerIdx = find(targIdx==0);
    Header{headerIdx, 1} = ['-SamplingFrequency ' num2str(newFs)];

    %% Create new file name:
    CSCFilename = ['SampleRate' num2str(newFs) 'Hz_' CSCFilename ];
    

end

%% Save the results to a new Neuralynx continuous file (.NCS):
if ~isequal(OriginalCSCFilename, CSCFilename)
    % Create necessary input vectors for the new file:
    ChannelNumbers = chanNum * ones(1,length(Timestamps));
    SampleFrequencies = newFs * ones(1,length(ChannelNumbers));
    NumberOfValidSamples = numValidSamp * ones(1,length(Timestamps));
    
    cscFile = fullfile(CSCFilePath, CSCFilename);
    Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
        SampleFrequencies, NumberOfValidSamples, Samples, Header);
end
clear all


%% Sub-function to fill in time gaps if recording was paused:
function [samples, timeStamps] = fillTimeGaps(gapIdx, samples, timeStamps)
numBreaks = size(gapIdx,2);
gapIdx = [0 gapIdx];
newTimeStamps = [];
newSamples = [];
nsamp = 512;

for j = 1:numBreaks
    tempTS = timeStamps((gapIdx(j)+1):gapIdx(j+1));
    tempSamp = samples(:,(gapIdx(j)+1):gapIdx(j+1));
    precise_samples=tempSamp(:);
    clear tempSamp
    eelen=length(tempTS);
    precise_timestamps = zeros(eelen*nsamp, 1);
    idx = 1;
    for i = 1:eelen
      if i < eelen
        t1 = tempTS(i);
        t2 = tempTS(i+1);
        interval = (t2-t1)/nsamp;
        trange =([t1 : interval : t2-interval]);
        precise_timestamps(idx:idx+nsamp-1,1) = trange;
      else
        t1 = tempTS(i);
        t2 = t1+interval*nsamp;
        trange =([t1 :interval : t2-interval]);
        precise_timestamps(idx:idx+nsamp-1,1) = trange;
      end
      idx = idx + nsamp;
    end
    gapTime = timeStamps(gapIdx(j+1)+1) - precise_timestamps(end);
    numFillTs = floor(gapTime/interval - 1);
    fillTime = interval:interval:(numFillTs*interval);
    fillTime = precise_timestamps(end) + fillTime;
    fillTime = fillTime(:);
    fillSamp = fillTime .* 0;
    newTimeStamps = [newTimeStamps; precise_timestamps; fillTime];
    newSamples = [newSamples; precise_samples; fillSamp];
end

tempTS = timeStamps((gapIdx(end)+1):end);
tempSamp = samples(:,(gapIdx(end)+1):end);
precise_samples=tempSamp(:);
clear tempSamp
eelen=length(tempTS);
precise_timestamps = zeros(eelen*nsamp, 1);
idx = 1;
for i = 1:eelen
  if i < eelen
    t1 = tempTS(i);
    t2 = tempTS(i+1);
    interval = (t2-t1)/nsamp;
    trange =([t1 : interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  else
    t1 = tempTS(i);
    t2 = t1+interval*nsamp;
    trange =([t1 :interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  end
  idx = idx + nsamp;
end

newTimeStamps = [newTimeStamps; precise_timestamps];
newSamples = [newSamples; precise_samples];
newEnd = floor(length(newSamples)/nsamp);
newLength = nsamp * newEnd;
timeStamps = newTimeStamps(1:nsamp:newLength)';
newSamples = newSamples(1:newLength);
samples = reshape(newSamples, nsamp, newEnd);

