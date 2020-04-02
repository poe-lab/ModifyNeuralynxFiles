function DetectSpikesInCSCs
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
mkdir(dataFolder, 'Detected_Spikes')
peakLocs = 8;
waveLength = 32;

for m = 1:numberOfDataFiles
    fileName = strtrim(fileList(m,:)); %Removes any white space at end of file name string.
    neuralynxFile = fullfile(dataFolder,fileName); %Full file path for Neuralynx file to be loaded
    newSpikeFile = fullfile(dataFolder, 'Detected_Spikes', ['spikes_' fileName]); %Full file path for new data file with reset time stamps
    newSpikeFile = strrep(newSpikeFile, '.ncs', '.nse');
    
    % Load the file to apply the wavelet filter to the signal:
    [Timestamps, ChannelNumbers, SampleFreq, numSamples, Samples, Header] =...
        Nlx2MatCSC(neuralynxFile, [1 1 1 1 1], 1, 1, [] );
    
    % Reshape data into a column vector:
    [nsamp,n1] = size(Samples);
    newM = nsamp*n1;
    Samples = reshape(Samples, newM, 1);
    
    %% Fill in the time stamps:
    eelen=length(Timestamps);
    if isempty(nsamp)
        nsamp=512;
    end
    precise_timestamps = zeros(eelen*nsamp, 1);
    idx = 1;
    for i = 1:eelen
      if i < eelen
        t1 = Timestamps(i);
        t2 = Timestamps(i+1);
        interval = (t2-t1)/nsamp;
        trange =([t1 : interval : t2-interval]);
        precise_timestamps(idx:idx+nsamp-1,1) = trange;
      else
        t1 = Timestamps(i);
        t2 = t1+interval*nsamp;
        trange =([t1 :interval : t2-interval]);
        precise_timestamps(idx:idx+nsamp-1,1) = trange;
      end
      idx = idx + nsamp;
    end 
    Timestamps = precise_timestamps;
    clear precise_timestamps
    
    %% Detect spikes
    % Zoran Nenadic spike detection method:
    %spikesIdx = detect_spikes_wavelet(Samples',SampleFreq/1000,[0.5 1.5],11,'l',0.1,'bior1.5',0,0);
    
    % Quian Quiroga spike detection method:
    noiseEstStdDev = median(abs(Samples)/0.6745);
    threshold = 5 * noiseEstStdDev; % I originally used 4x because that was
    % what was in Quiroga et al (Neural Computation, 2004). Change to 5x
    % which what Quiroga has in Wave_Clus v2.0 (2009).
    spikesIdx = [];

    spikesIdx = find(Samples > threshold);
    for j = 1:5
        e=diff(spikesIdx) > j-1;
        e=[true; e];
        spikesIdx = spikesIdx(e);
    end
    Y=[];
    a=spikesIdx-7;
    b=spikesIdx+7;
    for j = 1:size(spikesIdx,1)
        Y(j,:)=Samples(a(j):b(j));
    end
    Y = Y';
    [~,I] = max(Y);
    spikesIdx = spikesIdx + I' - 8;
    for j = 1:8
        e=diff(spikesIdx) > j-1;
        e=[true; e];
        spikesIdx = spikesIdx(e);
    end

    spikeTimesClear = Timestamps(spikesIdx);
    clear Timestamps

    [spikeData, spikeTs] = extractWaveforms05312017(Samples', spikesIdx, spikeTimesClear, peakLocs, waveLength);
    clear Samples qSpikesIdx spikeTimesClear
    % Find File UUID in Header:
    targ= strfind(Header,'-FileUUID');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    FileUuidIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Find Session UUID in Header:
    targ= strfind(Header,'-SessionUUID');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    SessionUuidIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Find time created in Header:
    targ= strfind(Header,'-TimeCreated');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    TimeCreatedIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Find time closed in Header:
    targ= strfind(Header,'-TimeClosed');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    TimeClosedIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Find reference channel in Header:
    targ= strfind(Header,'-Reference');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    ReferenceIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Find AD bit volts in Header:
    targ= strfind(Header,'-ADBitVolts');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    ADBitVoltsIdx = find(targIdx==0);   
    clear targ targIdx
    ADBitVolts = str2double(strrep(Header{ADBitVoltsIdx,1}, '-ADBitVolts ', ''));
    uVThreshold = threshold * ADBitVolts * 1000000; % convert threshold to microvolts
    % Find AD channel in Header:
    targ= strfind(Header,'-ADChannel');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    ADChannelIdx = find(targIdx==0);   
    chNum = str2double(strrep(Header{ADChannelIdx,1}, '-ADChannel ', ''));
    clear targ targIdx
    
    % Find input range in Header:
    targ= strfind(Header,'-InputRange');
    for i=1:length(targ)
        targIdx(i)= isempty(targ{i});
    end
    InputRangeIdx = find(targIdx==0);   
    clear targ targIdx
    
    % Customize spike header
    spikeHeader = {'######## Neuralynx Data File Header';
        '-FileType Spike';'-FileVersion 3.4';
        Header{FileUuidIdx,1};
        Header{SessionUuidIdx,1};
        '-ProbeName '; ['-OriginalFileName "' newSpikeFile '"'];
        Header{TimeCreatedIdx,1};Header{TimeClosedIdx,1}; '';
        '-RecordSize 112'; '-ApplicationName Cheetah "6.0.1 "';
        '-AcquisitionSystem AcqSystem1 DigitalLynxSX';
        Header{ReferenceIdx,1}; '-SamplingFrequency 32000';
        '-ADMaxValue 32767'; Header{ADBitVoltsIdx,1};
        ['-AcqEntName SE' num2str(chNum+1)]; '-NumADChannels 1';
        Header{ADChannelIdx,1};  Header{InputRangeIdx,1};
        '-InputInverted True'; '';
        '-DSPLowCutFilterEnabled False'; '-DspLowCutFrequency 300'; % Not accurate
        '-DspLowCutNumTaps 256'; '-DspLowCutFilterType FIR';
        '-DSPHighCutFilterEnabled False'; '-DspHighCutFrequency 6000';
        '-DspHighCutNumTaps 256'; '-DspHighCutFilterType FIR';
        '-DspDelayCompensation False'; '-DspFilterDelay_µs 3984';
        '';'-WaveformLength 32'; '-AlignmentPt 8';
        ['-ThreshVal ' num2str(round(uVThreshold))];
        '-MinRetriggerSamples 8'; '-SpikeRetriggerTime 250';
        '-DualThresholding False'; '';
        '-Feature Peak 0 0 0 31 1'; '-Feature Valley 1 0 0 31 1';
        '-Feature Energy 2 0 0 31 1'; '-Feature Height 3 0 0 31 1';
        '-Feature NthSample 4 0 0 31 1 4'; '-Feature NthSample 5 0 0 31 1 16';
        '-Feature NthSample 6 0 0 31 1 24'; '-Feature NthSample 7 0 0 31 1 28'};
    
    %% Reshape waveform data:
    spikeData = permute(spikeData,[2 3 1]);

    %% Create Features variable:
    X = min(spikeData, [],1);
    X = squeeze(X);
    Y= max(spikeData,[],1);
    Y = squeeze(Y);
    filler = zeros(size(Y,1),6);
    Q = [Y X filler]';
    Features = Q;
    clear X Y Q
    ScNumbers = chNum * ones(1, length(spikeTs));
    cellNumbers = zeros(1, length(spikeTs));
    Mat2NlxSpike(newSpikeFile, 0, 1, [], [1 1 1 1 1 1], spikeTs',...
        ScNumbers, cellNumbers, Features, spikeData, spikeHeader);
    clear spikeData spikeTs newSpikeFile ScNumbers cellNumbers Features spikeHeader
end