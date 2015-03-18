function OverlaySpikes(F, spikes,varargin)

SAVE = 0;
SAVEPATH = '';

% ------------------------------------------------------------------------------
% parse varargin
% ------------------------------------------------------------------------------
if nargin > 1
    for i = 1:2:length(varargin)-1
        if isnumeric(varargin{i+1})
            eval([varargin{i} '= [' num2str(varargin{i+1}) '];']);
        else
            eval([varargin{i} '=' char(39) varargin{i+1} char(39) ';']);
        end
    end
end

figure
t = 0:1/10:length(F)/10-1/10;
plot(t,F,'k');
Add_Event_Lines(gca,spikes/10);

if(SAVE)
    print('-djpeg',SAVEPATH);
end