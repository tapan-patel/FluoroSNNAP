function [Nhigh_activity_frames,ensemble_per_second,ensemble_frames,CoreEnsembles,EnsembleR,spk]=NetworkEnsemble(s)
% Compute network ensembles
% Input: data structure generated from PostProcess.m
% Reference: Visual stimuli recruit intrinsically generated cortical ensembles.
% Miller JE, Ayzenshtat I, Carrillo-Reid L, Yuste R.
% PNAS 111(38):E4053-61 2014

% Author: Tapan P Patel, PhD, tapan.p.patel@gmail.com
% University of Pennsylvania
load('params.mat');
multiWaitbar('Detecting network ensembles',0,'Name','');
if(~isfield(s,'foopsi') || isempty(s.foopsi))
    % Fast spike inference
    spk = zeros(size(s.dF_cell));
    
    for m=1:s.N
        
        x = run_oopsi(s.dF_cell(m,:));
        spk(m,:) = x.n';
        spk(m,:) = spk(m,:)./max(spk(m,:));
    end
    s.foopsi = spk;
end
spk = s.foopsi;

try
    
    thresh = params.network_ensemble_sd;
    Nsur = params.network_ensemble_Nsur;
catch
    errordlg('Parameters for network ensemble not set. Please try again.','Bad input','modal');
    return;
end
[M,T] = size(spk);
% APs is binary activity data
APs = zeros(size(spk));
for i=1:size(spk,1)
    sd = std(spk(i,:));
    APs(i,:) = spk(i,:)>sd*thresh;
end
Nactive_cells = zeros(Nsur,T); % Total number of active cells per frame. Each row is a separate surrogate resampling experiment
% Shuffle binary activity data 1000 times
for i=1:Nsur
    multiWaitbar('Detecting network ensembles',i/Nsur*.5);
    APsur = zeros(size(APs));
    for k=1:size(APs,1)
        APsur(k,:) = APs(k,randperm(length(APs(k,:))));
    end
    Nactive_cells(i,:) = sum(APsur);
end

% For each frame, compare the total number of actual active neurons to
% number of surrogate active neurons. If actual is >95 percentile of
% surrogate, this is p<0.05 significant for co-active neurons
coactive_thresh = zeros(1,T);
for i=1:T
    coactive_thresh(i) = prctile(Nactive_cells(:,i),95);
end
Ncoactive = sum(APs);
coactive_frames = Ncoactive>coactive_thresh;

Nhigh_activity_frames = sum(coactive_frames);
ensemble_per_second = Nhigh_activity_frames/T*s.fps;

%% Correlated ensembles
% The spatial similarity (i.e. groups of neurons that are co-activated)
% between ensembles is evaluated using Pearson's correlation coefficient
EnsembleR = zeros(Nhigh_activity_frames);
ensemble_frames = find(coactive_frames);
for i=1:Nhigh_activity_frames
    for j=i+1:Nhigh_activity_frames
        if(i~=j)
            EnsembleR(i,j) = corr(APs(:,ensemble_frames(i)),APs(:,ensemble_frames(j)));
            EnsembleR(j,i) = EnsembleR(i,j);
        end
    end
end
EnsembleR(isnan(EnsembleR))=0;

% Statistical testing: For each pair-wise ensemble,
% generate 1000 independent surrogate ensembles by randomizing active
% cells while preserving the total number of active cells per frame in one of the frames, i.e.
% shuffling across cells. Compute pair-wise correlation coefficient for
% each pair-wise surrogate ensemble. Threshold for significance = r > 95th
% percentile of surrogate dataset
P = zeros(Nhigh_activity_frames);
Nsur_ensemble = Nsur;
for i=1:Nhigh_activity_frames
    multiWaitbar('Detecting network ensembles',i/Nhigh_activity_frames*.5+.5);
    if(params.parallel)
        parfor j=1:Nhigh_activity_frames
            if(i~=j)
                Rsur = zeros(Nsur_ensemble,1);
                for k=1:Nsur_ensemble
                    % Generate 1000 surrogate ensembles by shuffling active cells in
                    % frame j. Compute pair-wise r between i and surrogate j
                    
                    x = APs(:,ensemble_frames(j));
                    sur = x(randperm(M));
                    Rsur(k) = corr(APs(:,ensemble_frames(i)),sur);
                end
                Rsur(isnan(Rsur))=0;
                thr = prctile(Rsur,95);
                if(EnsembleR(i,j)>thr)
                    P(i,j) = 1;
                    %             P(j,i) = 1;
                end
            end
        end
    else
        for j=1:Nhigh_activity_frames
            if(i~=j)
                Rsur = zeros(Nsur_ensemble,1);
                for k=1:Nsur_ensemble
                    % Generate 1000 surrogate ensembles by shuffling active cells in
                    % frame j. Compute pair-wise r between i and surrogate j
                    
                    x = APs(:,ensemble_frames(j));
                    sur = x(randperm(M));
                    Rsur(k) = corr(APs(:,ensemble_frames(i)),sur);
                end
                Rsur(isnan(Rsur))=0;
                thr = prctile(Rsur,95);
                if(EnsembleR(i,j)>thr)
                    P(i,j) = 1;
                    %             P(j,i) = 1;
                end
            end
        end
    end
end

%% Core-ensembles
% Defined as a group of coactive neurons that are conserved in all
% significantly correlated ensembles
CoreEnsembles = cell(0);

for i=1:Nhigh_activity_frames
    
    idx = find(P(i,:));
    if(~isempty(idx))
        % Find the intersect or common neurons between idx
        common_neurons = intersect(find(APs(:,ensemble_frames(i))),find(APs(:,ensemble_frames(idx(1)))));
        for j=2:length(idx)
            common_neurons = intersect(common_neurons,find(APs(:,ensemble_frames(idx(j)))));
        end
        if(~isempty(common_neurons))
            CoreEnsembles(end+1) = {common_neurons};
        end
    end
end

