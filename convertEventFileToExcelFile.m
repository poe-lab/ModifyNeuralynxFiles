%function convertEventFileToExcelFile
Filename = 'C:\Users\Brooks\UMICH_GoogleDrive\AnimalResearch\RatDataAnalysis\Rat1026\Events.nev';
FieldSelection = [1 0 0 0 1];
ExtractHeader = 0;
ExtractMode = 1;
ModeArray = [];
[TimeStamps, EventStrings] = Nlx2MatEV( Filename, FieldSelection, ExtractHeader, ExtractMode, ModeArray );



a=0;
b=1;
while isequal(a, 0)
    if (Aprime(b+1) - Aprime(b)) < 5
        Aprime(b+1) = [];
    else
        b = b+1;
    end
    c = size(Aprime,1);
    if b+1>c
        a =1
    end
        
end