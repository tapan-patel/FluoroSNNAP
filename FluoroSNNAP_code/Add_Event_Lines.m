% function [] = Add_Event_Lines(fh, Events, Event_Labels)
%
% Inputs:
%       fh  - figure handle of axes to add event lines and labels
%       Events  - t vector of time to draw lines
%       Event_Labels    - single entry or list of text labesl corresponding to events
%
% Outputs:
%
% Template borrowed from Anthony Choo. 
%
% ------------------------------------------------------------------------------
function [] = Add_Event_Lines(hax, Events,varargin)

% ------------------------------------------------------------------------------
% Default Arguments
% ------------------------------------------------------------------------------
C = 'b';
LW = 1;
FS = 13;

% ------------------------------------------------------------------------------
% parse varargin
% ------------------------------------------------------------------------------
if nargin > 2
    for i = 1:2:length(varargin)-1
        if isnumeric(varargin{i+1})
            eval([varargin{i} '= [' num2str(varargin{i+1}) '];']);
        else
            eval([varargin{i} '=' char(39) varargin{i+1} char(39) ';']);
        end
    end
end



Fig_ylim = get(hax,'ylim');


for j = 1:length(Events)
    lh=line([Events(j) Events(j)],Fig_ylim);
    set(lh,'linestyle','-','color',C,'linewidth',LW);
end

