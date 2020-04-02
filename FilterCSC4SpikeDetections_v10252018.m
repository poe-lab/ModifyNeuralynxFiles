function FilterCSC4SpikeDetections_v10252018
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
mkdir(dataFolder, 'FilterCSC4SpikeDetection')

for m = 1:numberOfDataFiles
    CSCFilename = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    cscFile = fullfile(dataFolder,CSCFilename); %Full file path for Neuralynx file to be loaded
    
    % load data:    
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
    % This IF statement catches Cheetah 160 recordings.
    if isempty(LowFreqCutIdx)
        clear targ targIdx
        targ= strfind(Header,'-AmpLowCut');
        for i=1:length(targ)
            targIdx(i)= isempty(targ{i});
        end
        LowFreqCutIdx = find(targIdx==0);
        recordHighPass = str2double(strrep(Header{LowFreqCutIdx,1}, '-AmpLowCut ', '')); %#ok<*FNDSB>
        cheetah160 = 1;
    else
        recordHighPass = str2double(strrep(Header{LowFreqCutIdx,1}, '-DspLowCutFrequency ', ''));
        cheetah160 = 0;
    end 
    clear targ targIdx
    
    %% Find low pass filter in Header:
    targ= strfind(Header,'-DspHighCutFrequency');

    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    HighFreqCutIdx = find(targIdx==0);   
    % This IF statement catches Cheetah 160 recordings.
    if isempty(HighFreqCutIdx)
        clear targ targIdx
        targ= strfind(Header,'-AmpHiCut');
        for i=1:length(targ)
            targIdx(i)= isempty(targ{i});
        end
        HighFreqCutIdx = find(targIdx==0);
        recordLowPass = str2double(strrep(Header{HighFreqCutIdx,1}, '-AmpHiCut ', ''));
    else
        recordLowPass = str2double(strrep(Header{HighFreqCutIdx,1}, '-DspHighCutFrequency ', ''));
    end 
    clear targ targIdx
    
    

    %% Design the filter:
    highPassFreq = 600;
    lowPassFreq = [];
    notchSetting = 0;
    newSampFreq = Fs;
    [highPassFreq, ~, sos, g, ~] = filterSettingsCheck(highPassFreq, lowPassFreq, recordHighPass, recordLowPass, Fs, newSampFreq, notchSetting);


    %% Filter the data:
    if ~isempty(g)
        Samples = filtfilt(sos,g, Samples);
        
        % Create new file name:
        CSCFilename = ['Filtered_' CSCFilename ];
        
        % Update High Cut and Low Cut frequencies in Header:
        if cheetah160 == 0 % Recorded on Digital Lynx system
%             Header{HighFreqCutIdx,1} = ['-DspHighCutFrequency ' num2str(lowPassFreq)];
            Header{LowFreqCutIdx,1} = ['-DspLowCutFrequency ' num2str(highPassFreq)];
        else % Recorded on the Cheetah 160 system
%             Header{HighFreqCutIdx,1} = ['-AmpHiCut ' num2str(lowPassFreq)];
            Header{LowFreqCutIdx,1} = ['-AmpLowCut ' num2str(highPassFreq)];
        end
    end

    %% Save the results to a new Neuralynx continuous file (.NCS):
    if ~isequal(OriginalCSCFilename, CSCFilename)
        % Reshape the signal to Neuralynx format:
        [m2,~]=size(Samples);
        newEnd = floor(m2/512);
        Samples = reshape(Samples, m1, newEnd);

        % Create necessary input vectors for the new file:
        ChannelNumbers = chanNum * ones(1,length(Timestamps));
        SampleFrequencies = newSampFreq * ones(1,length(ChannelNumbers));
        NumberOfValidSamples = numValidSamp * ones(1,length(Timestamps));

        cscFile = fullfile(dataFolder, CSCFilename);
        Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ChannelNumbers,...
            SampleFrequencies, NumberOfValidSamples, Samples, Header);
    end
%     clear all
end

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
