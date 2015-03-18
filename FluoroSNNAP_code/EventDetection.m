function spks = EventDetection(dF)
% Input: fluorescence trace
% Output: indices of the onset of calcium transients
% Load fps from params.mat
try
    load('params.mat');
catch
    errordlg('params.mat does not exist. Go to Analyze -> Preferences -> Revert to default or set preferences','File does not exist','modal');
    return;
end
% params.event_type = 1 (template-based), 2 = deconvolution based
if(params.event_type==1)
    load('spikes.mat');
    spks = SpikesLib(dF',spikes);
else
    x = run_oopsi(dF);
    inferred_spks = x.n;
    idx = find(inferred_spks>params.event_thresh*std(inferred_spks));
    d=[2;diff(idx)];
    idx(d==1) = [];
    spks = idx;
    
end
spks(spks==1)=[];