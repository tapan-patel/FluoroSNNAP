function [spks, Call] = SpikesLib(x,spikes,varargin)
xorig = x;
spks = [];
try
    load('params.mat');
catch
    errordlg('Could not load params.mat. Run Analysis -> Preferences -> Revert to default.','File not found','modal');
    return
end
thr = params.event_thresh;
fps = params.fps;
if nargin > 2
    for i = 1:2:length(varargin)-1
        if isnumeric(varargin{i+1})
            eval([varargin{i} '= [' num2str(varargin{i+1}) '];']);
        else
            eval([varargin{i} '=' char(39) varargin{i+1} char(39) ';']);
        end
    end
end

% Spike library is built using 10fps acquisiton speed. If input signal is
% faster or slower, need to interpolate or downsample to match the library
if(fps>10)
    x = decimate(x,floor(fps/10));
elseif(fps<10)
    x = interp(x,floor(10/fps));
end
if(max(x)<10) % deltaF/F0 being used
    height = params.event_amplitude;
else % Either raw trace or intensity above background being used
    height = params.event_amplitude*min(x);
end
Call = zeros(length(spikes),length(x));
if(params.parallel)
    parfor i=1:length(spikes)
        snippet = spikes{i};
        L = length(snippet);
        C = zeros(size(x));
        for j=1:length(x)-(L+1)
            x_snippet = x(j:j+L-1);
            if(range(x_snippet)>height)
                R = corrcoef(x_snippet,snippet);
                C(j) = R(1,2);
            end
        end
        %     [~,idx]=findpeaks(C);
        %     spks = [spks idx(C(idx)>thr)];
        Call(i,:) = C;
        
    end
else
    for i=1:length(spikes)
        snippet = spikes{i};
        L = length(snippet);
        C = zeros(size(x));
        for j=1:length(x)-(L+1)
            x_snippet = x(j:j+L-1);
            if(range(x_snippet)>height)
                R = corrcoef(x_snippet,snippet);
                C(j) = R(1,2);
            end
        end
        %     [~,idx]=findpeaks(C);
        %     spks = [spks idx(C(idx)>thr)];
        Call(i,:) = C;
        
    end
end
spks = find(max(Call)>thr);
try
    spks([2*fps diff(spks)]<=1.5*fps)=[];
    spks = spks*fps/10;
end
% OverlaySpikes(xorig,spks);