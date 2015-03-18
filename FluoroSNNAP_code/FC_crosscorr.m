function [A,C] = FC_crosscorr(s)
% Functional connectivity using cross-correlation

% Pairwise correlations between neurons calculated using the maximum normalized cross-correlation.
% Onset of calcium events are used to generate a fluorescence traces, Fi and Fj for neurons i and j.
% The maximum normalized cross-correlation between i and j at lags 0 to
% 500ms is recorded.
% For statistical testing of functional connection between i and j,
% surrogate onset times are generated for j by random permutation.
% Surrogate fluorescence trace for j is computed. Normalized
% cross-correlation between Fi and Fsur_j is computed. This is done 1000
% times. If original max cross-correlation is >99 percentile of surrogate,
% a functional connection between i and j exists
%
%%
load('params.mat');
[N,T] = size(s.dF_cell);
A = zeros(N);
C = zeros(N);
Nsur = params.FC.CC.Nsur;
maxlag = params.FC.CC.maxlag; % 500 ms
for i=1:N
    multiWaitbar('Functional connectivity: cross-correlation',i/N);
    if(params.parallel)
        parfor j=1:N
            if(i~=j && ~isempty(s.Spikes_cell{i}) && ~isempty(s.Spikes_cell{j}))
                
                Csur = zeros(Nsur,1);
                Fi = Surrogate_Fluorescence(s.Spikes_cell{i},T,s.fps);
                Fj = Surrogate_Fluorescence(s.Spikes_cell{j},T,s.fps);
                r = xcorr(Fi,Fj,ceil(maxlag*s.fps),'coeff');
                C(i,j) = max(r);
                for k=1:Nsur
                    Nspks = length(s.Spikes_cell{j});
                    sur_spks = sort(randsample(T,Nspks));
                    F_sur_j = Surrogate_Fluorescence(sur_spks,T,s.fps);
                    r = xcorr(Fi,F_sur_j,ceil(maxlag*s.fps),'coeff');
                    Csur(k) = max(r);
                    
                end
                if(C(i,j)>prctile(Csur,99))
                    A(i,j) = 1;
                end
            end
        end
    else
        for j=1:N
            if(i~=j && ~isempty(s.Spikes_cell{i}) && ~isempty(s.Spikes_cell{j}))
                
                Csur = zeros(Nsur,1);
                Fi = Surrogate_Fluorescence(s.Spikes_cell{i},T,s.fps);
                Fj = Surrogate_Fluorescence(s.Spikes_cell{j},T,s.fps);
                r = xcorr(Fi,Fj,ceil(maxlag*s.fps),'coeff');
                C(i,j) = max(r);
                for k=1:Nsur
                    Nspks = length(s.Spikes_cell{j});
                    sur_spks = sort(randsample(T,Nspks));
                    F_sur_j = Surrogate_Fluorescence(sur_spks,T,s.fps);
                    r = xcorr(Fi,F_sur_j,ceil(maxlag*s.fps),'coeff');
                    Csur(k) = max(r);
                    
                end
                if(C(i,j)>prctile(Csur,99))
                    A(i,j) = 1;
                end
            end
        end
    end
    
end