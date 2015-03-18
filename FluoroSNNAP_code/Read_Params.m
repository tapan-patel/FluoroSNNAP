% Read_Params - return parameters for data processing
%
function P = Read_Params(Param_Name)


SAVE_RESULTS = 1;
EXTRACT_FIRST_FRAME = 0;

if(ispc)
%     SOURCE_PATH = 'C:\Users\Tapan\Documents\Research\Fluo4\';
%     TARGET_PATH = 'C:\Users\Tapan\Documents\Research\Fluo4\';
    SOURCE_PATH = 'J:\Tapan\';
    TARGET_PATH = 'J:\Tapan\';
else
    
    SOURCE_PATH = '/media/Windows/Users/Tapan/Documents/Research/Fluo4/';
    TARGET_PATH = '/media/Windows/Users/Tapan/Documents/Research/Fluo4/';
end

% ------------------------------------------------------------
% Analysis Options
% ------------------------------------------------------------

ANALYSIS_TYPE = 'BOTH'; % Options are WHOLE, CELL, BOTH or ExtractFirstFrame.

eval(['P=' Param_Name ';']);
