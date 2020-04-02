function CreateContinuousTTFileFromCSCs
for i = 1:4
    %% Pick each CSC file that is part of the tetrode:
    [CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
            'Pick CSC file.'},['Select CSC File #' num2str(i) ' of the Tetrode']);
    cscFile = fullfile(CSCFilePath, CSCFilename);
    [Timestamps,Samples, Header] = Nlx2MatCSC(cscFile, [1 0 0 0 1], 1, 1, [] );
    [m1,n1] = size(Samples);
    newM = m1*n1;
    Samples = reshape(Samples, newM, 1);
    
    %% Find the AD bit value in the header to convert to Volts:
    HedrLngth = length(Header);
    ADstr = 'ADBitVolts';
    for iHedr = 1:HedrLngth
        iADBV = strfind(Header{iHedr},ADstr);
        if isempty(iADBV)
        else
           iADBVcell = iHedr; 
            break
        end
    end
    ADBVstr = Header{iADBVcell};
    clear ADstr Header HedrLngth iHedr iADBV iADBVcell
    spltADBV = textscan(ADBVstr, '%s %f');
    ADBitVal = spltADBV{1,2};
    
    %% Convert amplitude values to Volts.
    data.STRM.data{i} = Samples*ADBitVal;
    clear ADBitVal ADBVstr spltADBV Samples
end

data.STRM.samprate=32556;
[filename,pathname] = uiputfile('TR_file.xls','Save training file as:');
save(filename, 'data')
clear all
end