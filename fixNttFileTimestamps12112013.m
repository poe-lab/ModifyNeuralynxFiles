function fixNttFileTimestamps12112013
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


for i = 1:3
    [NttFilename, NttFilePath] = uigetfile({'*.ntt',...
            'Pick joined NTT files.'},'Select Joined Spike File to Fix Header');
    nttFile = fullfile(NttFilePath, NttFilename);

    [Timestamps, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(nttFile, [1 1 1 1 1], 1, 1, [] );
    fixedNttFilename = ['fixedTimestamps' NttFilename];
    nttFile = fullfile(NttFilePath, fixedNttFilename);
    Timestamps = Timestamps-Timestamps(1,1) +1;
    if isequal(i,1)
    else 
        Timestamps = Timestamps + a;
    end
    a=Timestamps(1,end)+600000000;
%     Mat2NlxSpike(nttFile, 0, 1, [], [1 1 1 1 1 1], Timestamps, ScNumbers,...
%         CellNumbers, Features, Samples, Header);
end
clear all