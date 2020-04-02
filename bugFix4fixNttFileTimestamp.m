function bugFix4fixNttFileTimestamp
%Applies a fix to NTT timestamps generated in original version of
%'fixNttFileTimestamps where the first time point was t=0.  This will
%result in Neuralynx's Spike Sort 3D program not being able to load the
%spikes.
[NttFilename, NttFilePath] = uigetfile({'*.ntt',...
'Pick joined NTT files.'},'Select Joined Spike File to Fix Header');
nttFile = fullfile(NttFilePath, NttFilename);
[Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFile, [1 1 1 1 1], 1, 1, [] );
Timestamps = Timestamps + 1;
fixedNttFilename = ['tPlus1microsec' NttFilename];
nttFile = fullfile(NttFilePath, fixedNttFilename);
Mat2NlxSpike(nttFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers,...
CellNumbers, Features, Samples, Header);

