function [waveforms, spikeTs] = extractWaveforms05312017(data, spikeIdx, spikeTs, peakLocs, waveLength)
% Inputs:
%     data -- filtered data
%     spikeIdx -- indices of waveform peaks
%     spikeTs -- time stamps of waveform peaks
%     peakLocs -- desired location of the peak within waveform
%     waveLength -- desired number of points to extract per waveform

numWires = size(data,1);
numSamples = size(data,2);
modStart = spikeIdx > peakLocs;
spikeIdx = spikeIdx(modStart);
spikeTs = spikeTs(modStart);
clear modStart
modStop = spikeIdx < (numSamples + peakLocs - waveLength);
spikeIdx = spikeIdx(modStop);
spikeTs = spikeTs(modStop);

numSpikes = length(spikeIdx);
waveforms = zeros(numSpikes, waveLength, numWires);

for i = 1:numSpikes
    waveStart = spikeIdx(i) - peakLocs + 1;
    waveEnd = spikeIdx(i) - peakLocs + waveLength;
    waveforms(i, :, :) = data(:, waveStart : waveEnd)';
end

end


