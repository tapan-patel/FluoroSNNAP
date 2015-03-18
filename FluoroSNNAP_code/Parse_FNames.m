
function P = Parse_FNames(FNames, Tag)

switch Tag
    
    case 'Folder'
        P = FNames{1,1};
        
    case 'Exclude'
        P = FNames{1,2};
    case 'FrameRange'
        P = FNames{1,3};
        
    case 'Prefix'
        P = FNames{1,4};
%         

    otherwise
        disp('*** Parse_FNames unrecognized Tag');
        P = [];
end
