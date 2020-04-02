function RemoveEventsFromNEV
% SELECT DATA FILES:
%Request user input to select event file:
working_dir = pwd;

[filename, pathname] = uigetfile('*.nev', 'Select a Cheetah event file'); %This waits for user input and limits selection to .nev files.
% Check for whether or not a file was selected
if isempty(filename) || isempty(pathname)
uiwait(errordlg('You need to select an event file. Please try again',...
'ERROR','modal'));
cd(working_dir);
else
NevFile= fullfile(pathname, filename);
end
%load event file
ExtractHeader = 1;  % 0 for no and 1 for yes
ExtractMode = 1;  %Extract all data points
[eventTimestamps, EventIDs, TTLs, Extras, EventStrings, Header] =Nlx2MatEV( NevFile, [1 1 1 1 1], ExtractHeader, ExtractMode, []);

% Remove all events that do not pertain to air puffing:
k = strfind(EventStrings, 'TTL');

numEvents=size(k,1);
for i = 1:numEvents
    if isempty(k{i,1})
    else
        eventTimeStamps(i)=-1;
    end
end
m = eventTimeStamps>=0;
EventStrings = EventStrings(m);
eventTimeStamps = eventTimeStamps(m);
TTLs = TTLs(m);
Extras = Extras(:,m);
EventIDs = EventIDs(m);
AppendToFileFlag = 0;
ExportMode = 1;
ExportModeVector = [];
FieldSelectionFlags = [1 1 1 1 1 1];
modNevFile= fullfile(pathname, [filename(1:end-4) 'Modified.nev']);
Mat2NlxEV(modNevFile, AppendToFileFlag, ExportMode, ExportModeVector,...
              FieldSelectionFlags, eventTimestamps, EventIDs, TTLs, Extras, EventStrings, Header);
clear all
end

