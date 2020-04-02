function joinOverlappingNtt
working_dir = pwd;
% % Load first pre-joined NTT file's header:
% [NttFilename, NttFilePath] = uigetfile({'*.ntt',...
%         'Pick  1st pre-joined NTT files.'},'Select 1st Pre-Joined Spike Sorted Data File');
% cd(working_dir);
% nttFile = fullfile(NttFilePath, NttFilename);
% 
% [nttHeader1] = Nlx2MatSpike(nttFile, [0 0 0 0 0], 1, 1, [] );
% 
% % Load last pre-joined NTT file's header:
% [NttFilename, NttFilePath] = uigetfile({'*.ntt',...
%         'Pick  lastst pre-joined NTT.'},'Select Last Pre-Joined Spike Sorted Data File');
% nttFile = fullfile(NttFilePath, NttFilename);
% 
% [nttHeader2] = Nlx2MatSpike(nttFile, [0 0 0 0 0], 1, 1, [] );
% 
% % Create header for joined NTT file:
% nttHeader1{4,1} =  nttHeader2{4,1};



[NttFilename, NttFilePath] = uigetfile({'*.ntt',...
        'Pick first NTT file '},'Select 1st file to combine');
nttFile1 = fullfile(NttFilePath, NttFilename);
[Timestamps1, ScNumbers1, CellNumbers1, Features1, Samples1, Header1] =...
    Nlx2MatSpike(nttFile1, [1 1 1 1 1], 1, 1, [] );
lastUnit = max(CellNumbers1);

[NttFilename, NttFilePath] = uigetfile({'*.ntt',...
        'Pick 2nd NTT file '},'Select 2nd file to combine');
nttFile2 = fullfile(NttFilePath, NttFilename);
[Timestamps2, ScNumbers2, CellNumbers2, Features2, Samples2, Header2] =...
    Nlx2MatSpike(nttFile2, [1 1 1 1 1], 1, 1, [] );

CellNumbers2 = CellNumbers2 + lastUnit;

Timestamps = [Timestamps1, Timestamps2];
ScNumbers = [ScNumbers1, ScNumbers2];
CellNumbers = [CellNumbers1, CellNumbers2];
Features = [Features1, Features2];
Samples = cat(3, Samples1, Samples2);

[Timestamps,IX] = sort(Timestamps);
ScNumbers = ScNumbers(IX);
CellNumbers = CellNumbers(IX);
Features = Features(:,IX);
Samples = Samples(:, :, IX);

combinedNttFilename = ['Combined' NttFilename];
nttFile = fullfile(NttFilePath, combinedNttFilename);

Mat2NlxSpike(nttFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers,...
    CellNumbers, Features, Samples, Header1);

clear all