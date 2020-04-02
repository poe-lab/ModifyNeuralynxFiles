function modifyNEV_10012017
% Removes all events except for target event.
%% Request user input to select event file:
working_dir = pwd;
[filename, pathname] = uigetfile('*.nev', 'Select a Cheetah event file'); %This waits for user input and limits selection to .nev files.
% Check for whether or not a file was selected
if isempty(filename) || isempty(pathname)
    uiwait(errordlg('You need to select an event file. Please try again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    NevFile= fullfile(pathname, filename);
end
%% Load event file
[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] =...
    Nlx2MatEV(NevFile, [1 1 1 1 1], 1, 1, [] );

%% Remove all events that do not match target event:
k = strcmp(EventStrings, 'TTL Output on PCI-DIO24_0 board 0 port 0 value (0x0001).');
EventStrings = EventStrings(k);
Timestamps = Timestamps(k);
TTLs = TTLs(k);
Extras = Extras(:,k);
EventIDs = EventIDs(k);

%% Save modified NEV file:
modNevFile = fullfile(pathname, ['targetEvent_' filename]); %Full file path for new data file with reset time stamps
Mat2NlxEV(modNevFile, 0, 1, [], [1 1 1 1 1],...
	Timestamps, EventIDs, TTLs, Extras, EventStrings, Header);
end

