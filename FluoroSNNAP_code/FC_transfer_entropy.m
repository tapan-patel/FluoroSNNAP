function [maxte,ci] = FC_transfer_entropy(s)
% Computes functional connectivity using transfer entropy method. Returns
% peak transfer entropy across lags, concidence index (weighted transfer
% entropy across lags) and a binary adjacency matrix obtained by resampling
% spike times.
load('params');
Spikes_cell = s.Spikes_cell;

multiWaitbar('Functional connectivity: Transfer entropy',0);
%%
[M,T] = size(s.F_cell);
asdf = Spikes_cell;
asdf(end+1) = {s.fps};
asdf(end+1) = {[M,T]};
[maxte, ci] = ASDFTE(asdf, params.FC.TE.lags);

% removing self connections
maxte = maxte - diag(diag(maxte));
ci = ci - diag(diag(ci));
maxte_sur = [];
ci_sur = [];

%%
A = zeros(M); % Adjacency matrix
Nsur = params.FC.TE.Nsur; % number of times to reshuffle spike train for statistical testing.
% To test whether i is functionally connected to j, destroy temporal
% information from j and compute transfer entropy between i and j. If
% original > 95th percentile of surrogate, a functional connection exists
for i=1:M
    multiWaitbar('Functional connectivity: Transfer entropy',i/M);
    if(params.parallel)
        parfor j=1:M
            peakTE_sur = zeros(Nsur,1);
            CI_sur = zeros(Nsur,1);
            if(i~=j)
                for k=1:Nsur
                    spks = [];
                    if(~isempty(Spikes_cell(i)) && ~isempty(Spikes_cell(j)))
                        spks{1} = Spikes_cell{i};
                        spks{2} = randsample(M,numel(Spikes_cell{j}))';
                        
                        spks(end+1) = {s.fps};
                        spks(end+1) = {[2,3e3]};
                        [maxte_sur, ci_sur] = ASDFTE(spks,1:30);
                    end
                    peakTE_sur(k) = maxte_sur(1,2);
                    CI_sur(k) = ci_sur(1,2);
                end
                CI_sur(isnan(CI_sur)) = [];
                peakTE_sur(isnan(peakTE_sur)) = [];
                if(maxte(i,j) > prctile(peakTE_sur,95) || ci(i,j) > prctile(CI_sur,95))
                    A(i,j) = 1;
                end
            end
        end
    else
        for j=1:M
            peakTE_sur = zeros(Nsur,1);
            CI_sur = zeros(Nsur,1);
            if(i~=j)
                for k=1:Nsur
                    spks = [];
                    if(~isempty(Spikes_cell(i)) && ~isempty(Spikes_cell(j)))
                        spks{1} = Spikes_cell{i};
                        spks{2} = randsample(M,numel(Spikes_cell{j}))';
                        
                        spks(end+1) = {s.fps};
                        spks(end+1) = {[2,3e3]};
                        [maxte_sur, ci_sur] = ASDFTE(spks,1:30);
                    end
                    peakTE_sur(k) = maxte_sur(1,2);
                    CI_sur(k) = ci_sur(1,2);
                end
                CI_sur(isnan(CI_sur)) = [];
                peakTE_sur(isnan(peakTE_sur)) = [];
                if(maxte(i,j) > prctile(peakTE_sur,95) || ci(i,j) > prctile(CI_sur,95))
                    A(i,j) = 1;
                end
            end
        end
    end
end

%%
maxte = maxte.*A;
ci = ci.*A;