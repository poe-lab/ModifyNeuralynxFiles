function downsampleClosestto1kHzCSC_v01222019
% 01182019: Now supports multiple file selection for batch processing. -BAG
% 01222019: Fix for accurate time stamp interpolation. This is necessary
% for over-sampled files that had buffering issues (mostly for old Cheetah
% system). -BAG

%% Select all of the .NCS CSC files to be down-sampled: 
[fileSet,CSCFilePath] = uigetfile('*.ncs','Select Neuralynx CSC file(s)','MultiSelect', 'on');
numFiles = size(fileSet,2); % # of files selected

if iscell(fileSet)
    
else
    numFiles = 1;
end

%% Run the down-sampling program for each file:
for n = 1:numFiles
    %% Load the Neuralynx CSC (.NCS) file:
    if iscell(fileSet)
        CSCFilename = fileSet{1,n};
    else
        CSCFilename = fileSet;
    end
    cscFile = fullfile(CSCFilePath, CSCFilename);
    [TimeStamps, ChannelNumbers, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(cscFile, [1 1 1 1 1], 1, 1, [] );
    OriginalCSCFilename = CSCFilename;

    %% Reduce variables to only necessary information
    chanNum = ChannelNumbers(1);
    clear ChannelNumbers
    Fs = SampleFrequencies(1);
    clear SampleFrequencies
    nsamp = NumberOfValidSamples(1);
    clear NumberOfValidSamples

    %% Find and remove and timstamps = 0:
    timeZeroIdx = find(TimeStamps < 0);
    
    %% Check for out of order time stamps due to memory buffering error:
    modeDiffTS = [];
    testOrderTS = 1;
    while testOrderTS == 1    
        outOfOrderTS = find(diff(TimeStamps) < 0, 1);
        if isempty(outOfOrderTS)
            testOrderTS = 0;
        else
            if isempty(modeDiffTS)
                modeDiffTS = mode(diff(TimeStamps));
            end
            jumpToIdx = find(TimeStamps >= (TimeStamps(outOfOrderTS) + modeDiffTS), 1);
            TimeStamps= [TimeStamps(1:outOfOrderTS), TimeStamps(jumpToIdx:end)];
            Samples = [Samples(:, 1:outOfOrderTS), Samples(:, jumpToIdx:end)];
        end
    end
    
    %% Convert from binned to vector of samples:
    Samples=double(Samples(:)');
    Samples = Samples';
    % Find median sampling rate: 
    medianSampRate = (nsamp)/median(diff(TimeStamps)) * 10^6;
    
    %% Need to first interpolate the time stamps:
    timeInterval = 10^6 * 1/medianSampRate;
    eelen=length(TimeStamps);
    precise_timestamps = zeros(eelen*nsamp, 1);
    idx = 1;
    for i = 1:eelen
        t1 = TimeStamps(i);
        % Calculate ideal next time stamp in bins of nsamp length:
        t2 = TimeStamps(i) +  timeInterval*nsamp;
        % Interpolate time stamp for each sample in the bin:
        trange =(t1 : timeInterval : t2-timeInterval);
        % Place in the pre-allocated time vector:
        precise_timestamps(idx:idx+nsamp-1,1) = trange;
        % Modify index for next time bin:
        idx = idx + nsamp;
    end
    clear TimeStamps eelen idx t1 t2 trange
    
    %% Set a threshold:
    threshDiffTS = 2 * timeInterval;

    %% Fill in breaks in the data:
    timeBreaks = find(diff(precise_timestamps) > threshDiffTS);
    if isempty(timeBreaks)
        TimeStamps = precise_timestamps;
        clear precise_timestamps
    else
        numBreaks = size(timeBreaks,1);
        timeBreaks = [0; timeBreaks]; %#ok<*AGROW>
        newTimeStamps = [];
        newSamples = [];

        for j = 1:numBreaks
            % Get timestamps and samples up to next time break:
            tempTS = precise_timestamps((timeBreaks(j)+1):timeBreaks(j+1));
            tempSamp = Samples((timeBreaks(j)+1):timeBreaks(j+1));

            % Calculate length of break in microseconds:
            gapTime = precise_timestamps(timeBreaks(j+1)+1) - tempTS(end);

            % Calculate the number of time points to fill in break:
            numFillTs = floor(gapTime/timeInterval - 1);

            % Calculate the time stamps to fill in the gap:
            fillTime = timeInterval:timeInterval:(numFillTs*timeInterval);
            fillTime = tempTS(end) + fillTime;
            fillTime = fillTime(:);
            sizeFillTime = length(fillTime);
            if sizeFillTime > Fs
                % Fill with zeros if longer than 1 second break:
                fillSamp = zeros(sizeFillTime,1);
            else
                % Fill with NaN if < 1 second break to fill gaps later:
                fillSamp = NaN(sizeFillTime,1);
            end
            newTimeStamps = [newTimeStamps; tempTS; fillTime];
            newSamples = [newSamples; tempSamp; fillSamp];
        end
        %% Append the remaining data after the last break:
        newTimeStamps = [newTimeStamps; precise_timestamps(timeBreaks(numBreaks) + 1:end)];
        newSamples = [newSamples; Samples(timeBreaks(numBreaks) + 1:end)];
        
        clear timeBreaks numBreaks threshDiffTS numBreaks tempTS tempSamp...
            gapTime numFillTs fillTime precise_timestamps
        
        %% Fill gaps in signal using autoregressive modeling:
        Samples = fillgaps(newSamples, 500, 1);
        clear newSamples
        TimeStamps = newTimeStamps;
        clear newTimeStamps
        CSCFilename = ['FillGaps_'  CSCFilename];
    end

    
    %% Find high pass filter in Header:
    targ= strfind(Header,'-DspLowCutFrequency');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
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

    %% Calculate new sampling frequency closest to 1 kHz:
    downSampFactor = floor(Fs/1000);
    newSampFreq = Fs/downSampFactor;

    % Design the filter:
    highPassFreq = [];
    lowPassFreq = [];
    notchSetting = 0;
    [highPassFreq, lowPassFreq, sos, g, ~] = filterSettingsCheck(highPassFreq, lowPassFreq, recordHighPass, recordLowPass, Fs, newSampFreq, notchSetting);



    if downSampFactor == 1
        % If the sampling rate is as close to 1 kHz as possible, just convert
        % samples back to 512 x n for Neuralynx:
        [m2,~]=size(Samples);
        newEnd = floor(m2/512);
        shortEnd = newEnd * 512;
        Samples = Samples(1:shortEnd);
        Samples = reshape(Samples, m1, newEnd);
    else          
        %% Filter data if the current high cut filter *3 is > the down-sampled Fs:
        if newSampFreq < recordLowPass*3
            Samples = filtfilt(sos,g, Samples);
            dataHasBeenFiltered = 1;
        else
            dataHasBeenFiltered = 0;
        end
        %% Down-sample the signal and time stamps closest to 1000 samples/second:
%         newM = length(Samples) + mod(length(Samples),downSampFactor);
        Samples = Samples(1:downSampFactor:end);
        [m2,~]=size(Samples);
        newEnd = floor(m2/nsamp);
        shortEnd = newEnd * nsamp;
        Samples = Samples(1:shortEnd);
        TimeStamps = TimeStamps(1:downSampFactor:end);
        TimeStamps = TimeStamps(1:nsamp:end);
        TimeStamps = TimeStamps(1:newEnd);

        %% Filter the data if not done prior to downsampling:
        if dataHasBeenFiltered == 0 && ~isempty(g)
            Samples = filtfilt(sos,g, Samples);
        end

        %% Reshape the signal to Neuralynx format:
        Samples = reshape(Samples, nsamp, newEnd);

        %% Update High Cut and Low Cut frequencies in Header:
        if cheetah160 == 0 % Recorded on Digital Lynx system
            Header{HighFreqCutIdx,1} = ['-DspHighCutFrequency ' num2str(lowPassFreq)];
            Header{LowFreqCutIdx,1} = ['-DspLowCutFrequency ' num2str(highPassFreq)];
        else % Recorded on the Cheetah 160 system
            Header{HighFreqCutIdx,1} = ['-AmpHiCut ' num2str(lowPassFreq)];
            Header{LowFreqCutIdx,1} = ['-AmpLowCut ' num2str(highPassFreq)];
        end

        %% Update sampling frequency in Header:
        targ= strfind(Header,'-SamplingFreq');
        for i=1:length(targ)
            targIdx(i)= isempty(targ{i});
        end
        headerIdx = find(targIdx==0);
        Header{headerIdx, 1} = ['-SamplingFrequency ' num2str(newSampFreq)];

        %% Create new file name:
        CSCFilename = ['SampleRate' num2str(newSampFreq) 'Hz_' CSCFilename ];


    end

    %% Save the results to a new Neuralynx continuous file (.NCS):
    if ~isequal(OriginalCSCFilename, CSCFilename)
        % Create necessary input vectors for the new file:
        ChannelNumbers = chanNum * ones(1,length(TimeStamps));
        SampleFrequencies = newSampFreq * ones(1,length(ChannelNumbers));
        NumberOfValidSamples = nsamp * ones(1,length(TimeStamps));

        cscFile = fullfile(CSCFilePath, CSCFilename);
        Mat2NlxCSC(cscFile, 0, 1, [], [1 1 1 1 1 1], TimeStamps', ChannelNumbers,...
            SampleFrequencies, NumberOfValidSamples, Samples, Header);
    end
    clear TimeStamps ChannelNumbers SampleFrequencies NumberOfValidSamples Samples Header
end