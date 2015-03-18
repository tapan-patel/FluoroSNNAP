function [A,P] = FC_phase(s)
% Define adjacency matrix using departure of pair-wise phase difference
% from resampled phase
params = [];
load('params.mat');
if(~isfield(s,'phase') || isempty(s.phase))
    phi = zeros(size(s.F_cell));
    if(params.parallel)
        parfor i=1:size(s.F_cell,1)
            phi(i,:) = GetPhaseSpikes(s.Spikes_cell{i},size(s.F_cell,2));
        end
    else
        for i=1:size(s.F_cell,1)
            phi(i,:) = GetPhaseSpikes(s.Spikes_cell{i},size(s.F_cell,2));
        end
    end
    s.phase = phi;
end

[M,T] = size(s.phase);

A = zeros(M);
P = zeros(M);
Nsur = params.FC.phase.Nsur; % Number of times to resample
for i=1:M
    multiWaitbar('Functional connectivity: phase',i/M);
    if(params.parallel)
        parfor j=1:M
            if(i~=j)
                psur = zeros(Nsur,1);
                phi1 = s.phase(i,:);
                phi2 = s.phase(j,:);
                phi12 = mod((phi1-phi2),2*pi);
                for k=1:Nsur
                    
                    psi2 = AAFTsur(phi2,1);
                    
                    psi12 = mod((phi1-psi2'),2*pi);
                    [~,p] = kstest2(hist(phi12,linspace(0,2*pi,100)),...
                        hist(psi12,linspace(0,2*pi,100)));
                    psur(k) = p;
                end
                
                P(i,j) = mean(psur);
                
                if(prctile(psur,95)<params.FC.phase.alpha)
                    A(i,j) = 1;
                end
            end
        end
    else
        for j=1:M
            if(i~=j)
                psur = zeros(Nsur,1);
                phi1 = s.phase(i,:);
                phi2 = s.phase(j,:);
                phi12 = mod((phi1-phi2),2*pi);
                for k=1:Nsur
                    
                    psi2 = AAFTsur(phi2,1);
                    
                    psi12 = mod((phi1-psi2'),2*pi);
                    [~,p] = kstest2(hist(phi12,linspace(0,2*pi,100)),...
                        hist(psi12,linspace(0,2*pi,100)));
                    psur(k) = p;
                end
                
                P(i,j) = mean(psur);
                
                if(prctile(psur,95)<params.FC.phase.alpha)
                    A(i,j) = 1;
                end
            end
        end
    end
end
