function processed_analysis = PostProcess(varargin)
warning off
% Run this function after AnalyzeData.m to process the raw data - namely
% to detect spikes.

if(nargin>0)
    SOURCE_PATH = [];
    TARGET_PATH = [];
    FNames = varargin{1};
    if(nargin==2)
        files = varargin{2};
    end
else
    % Get basic parameters
    SOURCE_PATH = Read_Params('SOURCE_PATH');
    TARGET_PATH = Read_Params('TARGET_PATH');
    FNames = Read_FNames;
end
if(ispc)
    slash = '\';
else
    slash = '/';
end
try
    load('params.mat');
catch
    errordlg('Could not load params.mat. Please run Analyze -> Preferences -> Revert to defaults','File not found','modal');
    return;
end

%%
spikes = [];
load('spikes.mat');
for i = 1:length(FNames) % folders
    Folder = FNames{i};
    if(iscell(Folder))
        Folder = Folder{1};
    end
    
    if(~ispc)
        SOURCE_PATH = strrep(SOURCE_PATH,'\','/');
        Folder = strrep(Folder,'\','/');
    end
    %     Make sure that analysis.mat is present. If not, need to run AnalyzeData.m
    AnalysisPath = [SOURCE_PATH Folder slash 'analysis.mat'];
    AnalysisPath = strrep(AnalysisPath,'\','/');
    if(~exist(AnalysisPath,'file'))
        error(['analysis.mat does not exist in search path ' AnalysisPath '. Please run AnalyzeData.m first']);
    end
    load (AnalysisPath);
    files = {analysis.filename};
    
    disp('Stiching files together');
    % To stich the files, first find out the unique filename, ie those not
    % containing -file002.tif, etc
    unique_idx = find(cellfun(@(x) isempty(x),strfind(files,'file')))';
    unique_files = files(unique_idx);
    %j indexes continuation files
    for j=setdiff(1:numel(files),unique_idx)
        % Figure out the index of the parent file
        tmpidx = strfind(files{j},'file');
        prefix = [files{j}(1:tmpidx-2) '.tif'];
        % Now find where this parent file is
        parent_idx = find(strcmpi(unique_files,prefix));
        % Next figure out which continuation number this is, eg file002 is
        % #2, file003 is #3
        continuation_num = str2double(strtok(files{j}(tmpidx:end),['file','.tif']));
        
        % Finally, put this back in unique_idx
        unique_idx(parent_idx,continuation_num) = j;
    end
    processed_analysis = [];
    % We are ready to merge files
    for j=1:size(unique_idx,1)
        
        
        F_cell = [];
        F_whole = [];
        files_stiched = {};
        fps = [];
        Frames = zeros(nnz(unique_idx(j,:)),1);
        for k=1:nnz(unique_idx(j,:))
            if(isfield(analysis,'F_cell'))
                F = analysis(unique_idx(j,k)).F_cell;
                Frames(k) = size(F,2);
                F_cell = [F_cell,F];
            end
            if(isfield(analysis,'F_whole'))
                F = analysis(unique_idx(j,k)).F_whole;
                Frames(k) = size(F,2);
                F_whole = [F_whole;F];
            end
            
            files_stiched{k} = files{unique_idx(j,k)};
        end
        
        % Put it in a struct
        if(isfield(analysis,'F_cell'))
            processed_analysis(j).F_cell = F_cell;
        end
        if(isfield(analysis,'F_whole'))
            processed_analysis(j).F_whole = F_whole;
        end
        if(isfield(analysis,'image'))
            processed_analysis(j).image = analysis(unique_idx(j,1)).image;
        end
        processed_analysis(j).filename = files{unique_idx(j,1)};
        processed_analysis(j).files_stiched = files_stiched;
        processed_analysis(j).Frames = size(F_cell,2);
        
        if(isfield(analysis,'fps'))
            
            fps = analysis(unique_idx(j,1)).fps;
        else
            fps = 10;
            msgbox('FPS not set. Assuming acquisition rate of 10fps');
        end
        processed_analysis(j).fps = fps;
    end
    
    
    if(isfield(processed_analysis,'F_cell'))
        
        for j=1:numel(processed_analysis)
            % Open the segmentation file
            
            processed_analysis(j).L = analysis(j).L;
            
            % Write raw intensities to a txt file
            cprintf('*blue','%s\n', ['Writing raw data to ASCII-comma delimited file for ' processed_analysis(j).filename]);
            out_fid = [processed_analysis(j).filename(1:end-4) '.csv'];
            dlmwrite(out_fid,processed_analysis(j).F_cell');
        end
        
    end
    
    
    % Now we need to detect spikes, get transient info, do SCA, functional
    % connectivity, output results to txt file and make a figure
    
    for j=1:numel(processed_analysis)
        
        % Create a multipanel waitbar
        processed_analysis(j).params = params;
        
        % Initialize waitbars
        multiWaitbar('CloseAll','Name','');
        [~,y] = fileparts(processed_analysis(j).filename);
        multiWaitbar('Converting raw fluorescence to deltaF/F',0,'Name',y);
        multiWaitbar('Detecting onset of calcium events',0,'Name',y);
        if(params.analyze.sca)
            multiWaitbar('Synchronization Cluster Analysis',0);
        end
        if(params.analyze.FC)
            if(params.FC.method_idx==1)
                multiWaitbar('Functional connectivity: cross-correlation',0,'Name',y);
            elseif(params.FC.method_idx==2)
                multiWaitbar('Functional connectivity: partial correlation',0,'Name',y);
            elseif(params.FC.method_idx==3)
                multiWaitbar('Functional connectivity: phase',0,'Name',y);
            elseif(params.FC.method_idx==4)
                multiWaitbar('Functional connectivity: Granger causality',0,'Name',y);
            elseif(params.FC.method_idx==5)
                multiWaitbar('Functional connectivity: Transfer entropy',0,'Name',y);
            elseif(params.FC.method_idx==6)
                multiWaitbar('Functional connectivity: cross-correlation',0,'Name',y);
                multiWaitbar('Functional connectivity: partial correlation',0,'Name',y);
                multiWaitbar('Functional connectivity: phase',0,'Name',y);
                multiWaitbar('Functional connectivity: Granger causality',0,'Name',y);
                multiWaitbar('Functional connectivity: Transfer entropy',0,'Name',y);
            end
        end
        if(params.analyze.controllability)
            multiWaitbar('Analyzing network controllability & driver nodes',0);
        end
        if(params.analyze.kinetics)
            multiWaitbar('Computing calcium event transient kinetics',0);
        end
        if(params.analyze.spike_probability)
            multiWaitbar('Inferring spikes from fluorescence',0);
        end
        if(params.analyze.ensembles)
            multiWaitbar('Detecting network ensembles',0);
        end
        if(params.analyze.figure)
            multiWaitbar('Saving results and making summary figure',0);
        end
        
        
        
        if(isfield(processed_analysis,'F_whole'))
            F = processed_analysis(j).F_whole;
            processed_analysis(j).Spikes_whole = EventDetection(F);
        end
        if(isfield(processed_analysis,'F_cell'))
            F = processed_analysis(j).F_cell;
            [N,frames] = size(F);
            
            %%%%%%%%%%% Determine deltaF/F by subtracting each value with the
            %mean of the lower 50% of previous 10-s values and dividing it
            % by the mean of the lower 50% of previous 10-s values.
            
            dF_cell = zeros(size(F));
            
            interval_idx = unique([1:params.fps*params.F0_time:frames frames]);
            % For the first 10-s of data, take the min value for F0
            F0 = min(F(:,1:interval_idx(2)),[],2);
            for k=1:interval_idx(2)
                dF_cell(:,k) = (F(:,k)-F0)./F0;
            end
            for it=interval_idx(2):frames
                if mod(it,frames)==1
                    multiWaitbar('Converting raw fluorescence to deltaF/F',it/frames);
                end
                % x = previous 10-s of points
                x = F(:,it-params.F0_time*params.fps:it);
                pctl = prctile(x,50,2);
                F0 = zeros(N,1);
                for n=1:N
                    F0(n) = mean(x(n,x(n,:)<pctl(n)));
                end
                dF_cell(:,it) = (F(:,it)-F0)./F0;
            end
            
            multiWaitbar('Converting raw fluorescence to deltaF/F',1,'Color','g');
            processed_analysis(j).dF_cell = dF_cell;
            processed_analysis(j).N = N;
            
            Spikes_cell = cell(N,1);
            
            
            fprintf('\tDetecting onset of calcium transients\n');
            for k = 1:N % Get the spikes for all cells
                s =[];
                multiWaitbar('Detecting onset of calcium events',k/N,'Name',y);
                s = EventDetection(processed_analysis(j).dF_cell(k,:));
                Spikes_cell(k) = {s};
                
            end
            processed_analysis(j).Spikes_cell = Spikes_cell;
            
            multiWaitbar('Detecting onset of calcium events',1,'Color','g');
            
            
            % See if the user wants to do SCA
            if(params.analyze.sca)
            multiWaitbar('Synchronization Cluster Analysis',0);
            fprintf('\tPerforming synchronization analysis\n');
            C = SCA(processed_analysis(j),'DISPLAY',0,'save_flag',0,'wb',0);
            SI = max(C.SI);
            processed_analysis(j).SynchroCluster = C;
            processed_analysis(j).phase = C.phase;
            processed_analysis(j).SI = SI;
            multiWaitbar('Synchronization Cluster Analysis',1,'Color','g');
            end
            
            % See if the user wants to infer functional connectivity
            if(params.analyze.FC)
            % Functional connectivity
            fprintf('\tEstimating functional connectivity\n');
            if(params.FC.method_idx==1)
                [A,C] = FC_crosscorr(processed_analysis(j));
                clu = clustering_coef_bu(A);
                
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                
                processed_analysis(j).FC.CC.A = A;
                processed_analysis(j).FC.CC.C = C;
                processed_analysis(j).FC.CC.clustering_coef = clu;
                processed_analysis(j).FC.CC.modularity_Ci = Ci;
                processed_analysis(j).FC.CC.modularity_Q = Q;
                multiWaitbar('Functional connectivity: cross-correlation',1,'Color','g');
                
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.CC.Controllability.Nd = Nd;
                processed_analysis(j).FC.CC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.CC.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
                
            elseif(params.FC.method_idx==2)
                multiWaitbar('Functional connectivity: partial correlation','Busy');
                [A,rho] = FC_partialcorr(processed_analysis(j).dF_cell);
                clu = clustering_coef_bu(A);
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                processed_analysis(j).FC.PC.A = A;
                processed_analysis(j).FC.PC.rho = rho;
                processed_analysis(j).FC.PC.clustering_coef = clu;
                processed_analysis(j).FC.PC.modularity_Ci = Ci;
                processed_analysis(j).FC.PC.modularity_Q = Q;
                multiWaitbar('Functional connectivity: partial correlation',1,'Color','g');
                
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.PC.Controllability.Nd = Nd;
                processed_analysis(j).FC.PC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.PC.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
                
            elseif(params.FC.method_idx==3)
                [A,P] = FC_phase(processed_analysis(j));
                clu = clustering_coef_bu(A);
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                processed_analysis(j).FC.phase.A=A;
                processed_analysis(j).FC.phase.Pvals = P;
                processed_analysis(j).FC.phase.clustering_coef = clu;
                processed_analysis(j).FC.phase.modularity_Ci = Ci;
                processed_analysis(j).FC.phase.modularity_Q = Q;
                multiWaitbar('Functional connectivity: phase',1,'Color','g');
                
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.phase.Controllability.Nd = Nd;
                processed_analysis(j).FC.phase.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.phase.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
                
            elseif(params.FC.method_idx==4)
                [A,P,F] = FC_granger(processed_analysis(j).dF_cell);
                clu = clustering_coef_bd(A);
                try
                    [Ci,Q] = modularity_louvain_dir(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                [~,~,d] = degrees_dir(A);
                processed_analysis(j).FC.GC.A = A;
                processed_analysis(j).FC.GC.Pvals = P;
                processed_analysis(j).FC.GC.F = F;
                processed_analysis(j).FC.GC.clustering_coef = clu;
                processed_analysis(j).FC.GC.modularity_Ci = Ci;
                processed_analysis(j).FC.GC.modularity_Q = Q;
                multiWaitbar('Functional connectivity: Granger causality',1,'Color','g');
                
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.GC.Controllability.Nd = Nd;
                processed_analysis(j).FC.GC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.GC.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
                
            elseif(params.FC.method_idx==5)
                [peakTE,CI] = FC_transfer_entropy(processed_analysis(j));
                A = peakTE;
                clu = clustering_coef_wd(peakTE);
                try
                    [Ci,Q] = modularity_louvain_dir(peakTE);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                [~,~,d] = degrees_dir(peakTE);
                processed_analysis(j).FC.TE.peakTE = peakTE;
                 processed_analysis(j).FC.TE.A = peakTE;
                processed_analysis(j).FC.TE.CI = CI;
                processed_analysis(j).FC.TE.clustering_coef = clu;
                processed_analysis(j).FC.TE.modularity_Ci = Ci;
                processed_analysis(j).FC.TE.modularity_Q = Q;
                multiWaitbar('Functional connectivity: Transfer entropy',1,'Color','g');
                
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.TE.Controllability.Nd = Nd;
                processed_analysis(j).FC.TE.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.TE.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
                
            elseif(params.FC.method_idx==6)
                % Do all
                [A,C] = FC_crosscorr(processed_analysis(j));
                clu = clustering_coef_bu(A);
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                processed_analysis(j).FC.CC.modularity_Ci = Ci;
                processed_analysis(j).FC.CC.modularity_Q = Q;
                processed_analysis(j).FC.CC.A = A;
                processed_analysis(j).FC.CC.C = C;
                processed_analysis(j).FC.CC.clustering_coef = clu;
                multiWaitbar('Functional connectivity: cross-correlation',1,'Color','g');
                multiWaitbar('Analyzing network controllability & driver nodes','Busy');
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.CC.Controllability.Nd = Nd;
                processed_analysis(j).FC.CC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.CC.Controllability.Nconfigs = Nconfigs;
                
                
                multiWaitbar('Functional connectivity: partial correlation','Busy');
                [A,rho] = FC_partialcorr(processed_analysis(j).dF_cell);
                clu = clustering_coef_bu(A);
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                processed_analysis(j).FC.PC.modularity_Ci = Ci;
                processed_analysis(j).FC.PC.modularity_Q = Q;
                processed_analysis(j).FC.PC.clustering_coef = clu;
                processed_analysis(j).FC.PC.A = A;
                processed_analysis(j).FC.PC.rho = rho;
                multiWaitbar('Functional connectivity: partial correlation',1,'Color','g');
                
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.PC.Controllability.Nd = Nd;
                processed_analysis(j).FC.PC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.PC.Controllability.Nconfigs = Nconfigs;
                
                [A,P] = FC_phase(processed_analysis(j));
                clu = clustering_coef_bu(A);
                try
                    [Ci,Q] = modularity_louvain_und(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                d = sum(A);
                processed_analysis(j).FC.phase.modularity_Ci = Ci;
                processed_analysis(j).FC.phase.modularity_Q = Q;
                processed_analysis(j).FC.phase.A=A;
                processed_analysis(j).FC.phase.Pvals = P;
                processed_analysis(j).FC.phase.clustering_coef = clu;
                multiWaitbar('Functional connectivity: phase',1,'Color','g');
                
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.phase.Controllability.Nd = Nd;
                processed_analysis(j).FC.phase.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.phase.Controllability.Nconfigs = Nconfigs;
                
                [A,P,F] = FC_granger(processed_analysis(j).dF_cell);
                clu = clustering_coef_bd(A);
                try
                    [Ci,Q] = modularity_louvain_dir(A);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                [~,~,d] = degrees_dir(A);
                processed_analysis(j).FC.GC.modularity_Ci = Ci;
                processed_analysis(j).FC.GC.modularity_Q = Q;
                processed_analysis(j).FC.GC.clustering_coef = clu;
                processed_analysis(j).FC.GC.A = A;
                processed_analysis(j).FC.GC.Pvals = P;
                processed_analysis(j).FC.GC.F = F;
                multiWaitbar('Functional connectivity: Granger causality',1,'Color','g');
                
                [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
                processed_analysis(j).FC.GC.Controllability.Nd = Nd;
                processed_analysis(j).FC.GC.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.GC.Controllability.Nconfigs = Nconfigs;
                
                [peakTE,CI] = FC_transfer_entropy(processed_analysis(j));
                clu = clustering_coef_wd(peakTE);
                try
                    [Ci,Q] = modularity_louvain_dir(peakTE);
                catch
                    Ci = zeros(N,1);
                    Q = 0;
                end
                [~,~,d] = degrees_dir(peakTE);
                processed_analysis(j).FC.TE.modularity_Ci = Ci;
                processed_analysis(j).FC.TE.modularity_Q = Q;
                processed_analysis(j).FC.TE.peakTE = peakTE;
                processed_analysis(j).FC.TE.A = peakTE;
                processed_analysis(j).FC.TE.CI = CI;
                processed_analysis(j).FC.TE.clustering_coef = clu;
                multiWaitbar('Functional connectivity: Transfer entropy',1,'Color','g');
                
                [Nd,drivernodes,Nconfigs] = ExactControllability(peakTE,'plotting',0);
                processed_analysis(j).FC.TE.Controllability.Nd = Nd;
                processed_analysis(j).FC.TE.Controllability.drivernodes = drivernodes;
                processed_analysis(j).FC.TE.Controllability.Nconfigs = Nconfigs;
                multiWaitbar('Analyzing network controllability & driver nodes',1,'Color','g');
            end
            if(params.FC.method_idx==6) % Default to phase method if the user selects "all" for the purpose of making figure
                A = processed_analysis(j).FC.phase.A;
                clu = processed_analysis(j).FC.phase.clustering_coef;
            end
            
            
            % Sychronization by modules
            SI_m = zeros(max(Ci),1);
            for k=1:max(Ci)
                tmp_s.F_cell = processed_analysis(j).F_cell(Ci==k,:);
                tmp_s.Spikes_cell = processed_analysis(j).Spikes_cell(Ci==k);
                tmp_s.fps = processed_analysis(j).fps;
                c_mod = SCA(tmp_s);
                SI_m(k) = max(c_mod.SI);
            end
            processed_analysis(j).SI_m = SI_m;
             processed_analysis(j).modules = Ci;
            processed_analysis(j).modularity = Q;
            end
            % Mean oscillation period (s)
            ISI = cell(N,1);
            for k=1:N
                ISI{k} = diff(processed_analysis(j).Spikes_cell{k});
            end
            
            OP = zeros(N,1);
            for k=1:N
                try
                    OP(k) = mean(ISI{k})/processed_analysis(j).fps;
                end
            end
              processed_analysis(j).OP = OP;
            % See if the user wants to determine transient kinetics
            if(params.analyze.kinetics)
            fprintf('\tDetermining transient kinetics - amplitude, rise time, and fall time\n');
            
            [DF,rise_time,fall_time,CV] = GetTransients(processed_analysis(j));
            multiWaitbar('Computing calcium event transient kinetics',1,'Color','g');
            ci = bootci(40,@(x) mean(x),DF(~isnan(DF)));
            
            Rise = [rise_time{:}]; Fall = [fall_time{:}];
            ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));
            
            ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
             processed_analysis(j).DF = DF;
             processed_analysis(j).rise_time = rise_time;
            processed_analysis(j).fall_time = fall_time;
            processed_analysis(j).CV = CV;
            end
            
            % See if the user wants to infer spike probability
            if(params.analyze.spike_probability)
            fprintf('\tPerforming spike inference\n');
            % Fast spike inference
            spk = zeros(size(processed_analysis(j).F_cell));
            
            for m=1:processed_analysis(j).N
                multiWaitbar('Inferring spikes from fluorescence',m/processed_analysis(j).N);
                x = run_oopsi(processed_analysis(j).dF_cell(m,:));
                spk(m,:) = x.n';
                spk(m,:) = spk(m,:)./max(spk(m,:));
            end
            multiWaitbar('Inferring spikes from fluorescence',1,'Color','g');
              processed_analysis(j).foopsi = spk;
            end
          
            
            
            % Network ensembles
            if(params.analyze.ensembles)
                
            [Nhigh_activity_frames,ensemble_per_second,ensemble_frames,CoreEnsembles,EnsembleR,foopsi]=NetworkEnsemble(processed_analysis(j));
            multiWaitbar('Detecting network ensembles',1,'Color','g');
            processed_analysis(j).NetworkEnsemble.ensembles_per_second = ensemble_per_second;
            processed_analysis(j).NetworkEnsemble.ensemble_frames = ensemble_frames;
            processed_analysis(j).NetworkEnsemble.Nensembles = Nhigh_activity_frames;
            processed_analysis(j).NetworkEnsemble.CoreEnsembles = CoreEnsembles;
            processed_analysis(j).NetworkEnsemble.CorrelatedEnsembles = EnsembleR;
            if(~params.analyze.spike_probability)
                processed_analysis(j).foopsi = foopsi;
            end
            end
            
            
            % Create per neuron report in an table/
            % Baseline fluorescence, Number of events, <ISI>, <DF>, CV,rise time, fall time, CI, PI, # assemblies,
            % clustering coef
            dat = zeros(N,11);
            if(params.analyze.sca)
            C = processed_analysis(j).SynchroCluster;
            end
            for k=1:N
                x = sort(processed_analysis(j).F_cell(k,:));
                dat(k,1) = mean(x(1:ceil(.1*length(x))));
                dat(k,2) = numel(processed_analysis(j).Spikes_cell{k});
                dat(k,3) = OP(k);
                if(params.analyze.kinetics)
                dat(k,4) = DF(k);
                dat(k,5) = CV(k);
                r = rise_time{k}; tau = fall_time{k};
                dat(k,6) = mean(r(~isnan(r)));
                dat(k,7) = mean(tau(~isnan(tau)));
                end
                if(params.analyze.FC)
                dat(k,8) = sum(A(k,:))/N;
                end
                if(params.analyze.sca)
                dat(k,9) = max(C.PI(k,:));
                dat(k,10) = nnz(C.PI(k,:)>.01);
                end
                if(params.analyze.FC)
                dat(k,11) = clu(k);
                end
            end
            
            % Export everything to .txt file
            fid = fopen([processed_analysis(j).filename(1:end-4) '_summary.txt'],'w');
            if(fid)
                fprintf(fid,'Summary for file:\t %s\n\n',processed_analysis(j).filename);
                fprintf(fid,'Total cells (ROIs):\t %d\n',N);
                fprintf(fid,'Frames:\t %d\nframe rate:\t %d\nduration (s):\t %.02f\n',...
                    size(processed_analysis(j).F_cell,2),processed_analysis(j).fps,...
                    size(processed_analysis(j).F_cell,2)/processed_analysis(j).fps);
                % Output functional connectivity summaries
                if(params.FC.method_idx==1)
                    if(params.analyze.FC)
                    fprintf(fid,'Functional connectivity method: cross-correlation\n');
                    fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.CC.A))/N);
                    
                    fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.CC.modularity_Q);
                    fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.CC.modularity_Ci));
                    end
                    if(params.analyze.controllability)
                    fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.CC.Controllability.Nd);
                    fprintf(fid,'\tList of driver nodes:\t');
                    for q=1:length(processed_analysis(j).FC.CC.Controllability.drivernodes)
                        fprintf(fid,'%d, ',processed_analysis(j).FC.CC.Controllability.drivernodes(q));
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.CC.Controllability.Nconfigs);
                    end
                    end
                
                if(params.FC.method_idx==2)
                    if(params.analyze.FC)
                    fprintf(fid,'Functional connectivity method: partial-correlation\n');
                    fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.PC.A))/N);
                    
                    fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.PC.modularity_Q);
                    fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.PC.modularity_Ci));
                    end
                    if(params.analyze.controllability)
                    fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.PC.Controllability.Nd);
                    fprintf(fid,'\tList of driver nodes:\t');
                    for q=1:length(processed_analysis(j).FC.PC.Controllability.drivernodes)
                        fprintf(fid,'%d, ',processed_analysis(j).FC.PC.Controllability.drivernodes(q));
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.PC.Controllability.Nconfigs);
                    end
                    end
                if(params.FC.method_idx==3)
                    if(params.analyze.FC)
                    fprintf(fid,'Functional connectivity method: phase\n');
                    fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.phase.A))/N);
                    
                    fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.phase.modularity_Q);
                    fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.phase.modularity_Ci));
                    end
                    if(params.analyze.controllability)
                    fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.phase.Controllability.Nd);
                    fprintf(fid,'\tList of driver nodes:\t');
                    for q=1:length(processed_analysis(j).FC.phase.Controllability.drivernodes)
                        fprintf(fid,'%d, ',processed_analysis(j).FC.phase.Controllability.drivernodes(q));
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.phase.Controllability.Nconfigs);
                    end
                    end
                
                if(params.FC.method_idx==4)
                    if(params.analyze.FC)
                    fprintf(fid,'Functional connectivity method: Granger causality\n');
                    fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.GC.A))/N);
                    
                    fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.GC.modularity_Q);
                    fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.GC.modularity_Ci));
                    end
                    if(params.analyze.controllability)
                    fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.GC.Controllability.Nd);
                    fprintf(fid,'\tList of driver nodes:\t');
                    for q=1:length(processed_analysis(j).FC.GC.Controllability.drivernodes)
                        fprintf(fid,'%d, ',processed_analysis(j).FC.GC.Controllability.drivernodes(q));
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.GC.Controllability.Nconfigs);
                    end
                    end
                if(params.FC.method_idx==5)
                    if(params.analyze.FC)
                    fprintf(fid,'Functional connectivity method: Transfer entropy\n');
                    fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.TE.A))/N);
                    
                    fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.TE.modularity_Q);
                    fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.TE.modularity_Ci));
                    end
                    if(params.analyze.controllability)
                    fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.TE.Controllability.Nd);
                    fprintf(fid,'\tList of driver nodes:\t');
                    for q=1:length(processed_analysis(j).FC.TE.Controllability.drivernodes)
                        fprintf(fid,'%d, ',processed_analysis(j).FC.TE.Controllability.drivernodes(q));
                    end
                    fprintf(fid,'\n');
                    fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.TE.Controllability.Nconfigs);
                    end
                end
                if(params.FC.method_idx==6)
                    if(params.analyze.FC)
                        fprintf(fid,'Functional connectivity method: cross-correlation\n');
                        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.CC.A))/N);
                        
                        fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.CC.modularity_Q);
                        fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.CC.modularity_Ci));
                        
                        fprintf(fid,'Functional connectivity method: partial-correlation\n');
                        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.PC.A))/N);
                        
                        fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.PC.modularity_Q);
                        fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.PC.modularity_Ci));
                        fprintf(fid,'Functional connectivity method: phase\n');
                        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.phase.A))/N);
                        
                        fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.phase.modularity_Q);
                        fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.phase.modularity_Ci));
                        fprintf(fid,'Functional connectivity method: Granger causality\n');
                        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.GC.A))/N);
                        
                        fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.GC.modularity_Q);
                        fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.GC.modularity_Ci));
                        fprintf(fid,'Functional connectivity method: Transfer entropy\n');
                        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(processed_analysis(j).FC.TE.peakTE))/N);
                        
                        fprintf(fid,'\tModularity:\t %.02f\n',processed_analysis(j).FC.TE.modularity_Q);
                        fprintf(fid,'\tNumber of modules:\t %d\n',max(processed_analysis(j).FC.TE.modularity_Ci));
                        
                    end
                    if(params.analyze.controllability)
                        fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.CC.Controllability.Nd);
                        fprintf(fid,'\tList of driver nodes:\t');
                        for q=1:length(processed_analysis(j).FC.CC.Controllability.drivernodes)
                            fprintf(fid,'%d, ',processed_analysis(j).FC.CC.Controllability.drivernodes(q));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.CC.Controllability.Nconfigs);
                        
                        fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.PC.Controllability.Nd);
                        fprintf(fid,'\tList of driver nodes:\t');
                        for q=1:length(processed_analysis(j).FC.PC.Controllability.drivernodes)
                            fprintf(fid,'%d, ',processed_analysis(j).FC.PC.Controllability.drivernodes(q));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.PC.Controllability.Nconfigs);
                        fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.phase.Controllability.Nd);
                        fprintf(fid,'\tList of driver nodes:\t');
                        for q=1:length(processed_analysis(j).FC.phase.Controllability.drivernodes)
                            fprintf(fid,'%d, ',processed_analysis(j).FC.phase.Controllability.drivernodes(q));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.phase.Controllability.Nconfigs);
                        fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.GC.Controllability.Nd);
                        fprintf(fid,'\tList of driver nodes:\t');
                        for q=1:length(processed_analysis(j).FC.GC.Controllability.drivernodes)
                            fprintf(fid,'%d, ',processed_analysis(j).FC.GC.Controllability.drivernodes(q));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.GC.Controllability.Nconfigs);
                        fprintf(fid,'\tNumber of driver nodes:\t %d\n',processed_analysis(j).FC.TE.Controllability.Nd);
                        fprintf(fid,'\tList of driver nodes:\t');
                        for q=1:length(processed_analysis(j).FC.TE.Controllability.drivernodes)
                            fprintf(fid,'%d, ',processed_analysis(j).FC.TE.Controllability.drivernodes(q));
                        end
                        fprintf(fid,'\n');
                        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',processed_analysis(j).FC.TE.Controllability.Nconfigs);
                        
                    end
                end
                
                if(params.analyze.sca)
                fprintf(fid,'Global synchronization index:\t %.02f\n',SI);
                end
                if(params.analyze.FC)
                fprintf(fid,'Synchronization index by modules:\t ');
                for k=1:length(SI_m)
                    fprintf(fid,'%f\t',SI_m(k));
                end
                fprintf(fid,'\n');
                end
                
                if(params.analyze.ensembles)
                
                fprintf(fid,'Number of network ensembles:\t %d\n',processed_analysis(j).NetworkEnsemble.Nensembles);
                fprintf(fid,'Ensembles per second:\t %.04f\n',processed_analysis(j).NetworkEnsemble.ensembles_per_second);
                for q=1:length(processed_analysis(j).NetworkEnsemble.CoreEnsembles)
                    fprintf(fid,'Core ensemble %d:\t',q);
                    ens = processed_analysis(j).NetworkEnsemble.CoreEnsembles{q};
                    for w=1:length(ens)
                        fprintf(fid,'%d, ',ens(w));
                    end
                    fprintf(fid,'\n');
                end
                end
                fprintf(fid,'Mean oscillation period (s):\t %.02f\n',mean(OP));
                
                if(params.analyze.kinetics)
                fprintf(fid,'Mean amplitude (deltaF/F0):\t %.04f\n95%% CI:\t%.04f\t %.04f\n',mean(DF(~isnan(DF))),ci(1),ci(2));
                fprintf(fid,'Coefficient of variation in amplitude:\t%f\n',mean(CV(~isnan(CV))));
                fprintf(fid,'Mean rise time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2));
                fprintf(fid,'Mean fall time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2));
                %Now output cell based report
                fprintf(fid,'\n\nNeuron ID\tBaseline fluorescence\tTotal events\t<ISI> (s)\tAmplitude\tCV\t<Rise time> (s)\t<Fall time> (s)\tConnectivity\tParticipation\t# assemblies\tClustering coef\n');
                for k=1:processed_analysis(j).N
                    fprintf(fid,'%d\t%f\t%.02f\t%.04f\t%.02f\t%f\t%.02f\t%.04f\t%.04f\t%d\t%.04f\t%0.04f\n',k,dat(k,1),dat(k,2),dat(k,3),dat(k,4),dat(k,5),dat(k,6),dat(k,7),dat(k,8),dat(k,9),dat(k,10),dat(k,11));
                end
                end
                fclose(fid);
                
                cprintf('*blue','%s\n', ['Summary statistics written to ' processed_analysis(j).filename(1:end-4) '_summary.txt']);
            else
                cprintf('*red','%s\n','ERROR: Could not write summary statistics. Unable to create file.');
            end
            if(params.analyze.figure)
            % Finally, make a composite figure and save it as eps
            multiWaitbar('Saving results and making summary figure','Busy');
            netfig = figure;
            set (netfig, 'Units', 'normalized', 'Position', [0,0,1,1],'PaperOrientation','portrait','PaperType','arch-C');
            subplot(2,2,1);
            if(isfield(processed_analysis(j),'image'))
                if(size(processed_analysis(j).image,3)~=3)
                    image = zeros([size(processed_analysis(j).image) 3]);
                    image(:,:,2) = imadjust(processed_analysis(j).image);
                elseif(size(processed_analysis(j).image,3)==3)
                    image = processed_analysis(j).image;
                else
                    image = logical(processed_analysis(j).L);
                end
            end
            imagesc(image);
            hold on
            
            c = regionprops(processed_analysis(j).L,'Centroid');
            c = reshape([c.Centroid],2,processed_analysis(j).N);
            x = c(1,:);
            y = c(2,:);
            Colors = varycolor(max(Ci));
            if(isempty(Colors))
                Colors = [0 0 1];
            end
            gplot(A,[x' y']);
            set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
            
            plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
            set(gca,'YDir','reverse');
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
            title({['N = ' num2str(processed_analysis(j).N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});
            
            subplot(2,2,2);
            [~,pos] = max(C.PI,[],2);
            [~,IDX] = sort(pos);
            imagesc(C.C(IDX,IDX),[0 1]); axis square; colorbar;
            colormap('jet');
            xlabel('Neuron ID (rearranged)'); ylabel('Neuron ID (rearranged)');
            msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
            title(msg);
            
            subplot(2,2,3);
            
            set(gca, 'ColorOrder', Colors);
            hold all;
            if(sum(Ci)==0)
                Ci = ones(N,1);
                scatter(x,y,'filled','SizeData',10^2);
            else
                for i=1:max(Ci)
                    scatter(x(Ci==i),y(Ci==i),'filled','SizeData',10^2);
                end
            end
            set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
            set(gca,'YDir','reverse');
            
            title_string = sprintf('N = %d. Modules = %d. Modularity = %f\n',length(x),max(Ci),Q);
            title_string = sprintf('%sSynchronization by modules: ',title_string);
            for k=1:length(SI_m)
                title_string = sprintf('%s %f,',title_string,SI_m(k));
            end
            title(title_string);
            
            subplot(2,2,4);
            t = 0:1/processed_analysis(j).fps:size(processed_analysis(j).F_cell,2)/processed_analysis(j).fps-1/processed_analysis(j).fps;
            hold all;
            cntr = 1;
            for q=1:max(Ci)
                modules = find(Ci==q);
                for m=1:length(modules)
                    try
                        evnts = processed_analysis(j).Spikes_cell{modules(m)};
                        plot(evnts/processed_analysis(j).fps,cntr,'.','Color',Colors(q,:));
                        cntr = cntr+1;
                    end
                end
            end
            title('Activity rearranged by modules');
            ylim([0 N]); xlim([0 frames/processed_analysis(j).fps]);
            xlabel('Time (s)'); ylabel('Neuron ID (rearranged)');
            % Save this figure in a new folder
            [pathstr, name] = fileparts(processed_analysis(j).filename);
            if(~exist([pathstr slash 'Figures'],'dir'))
                mkdir(pathstr,'Figures');
            end
            set(netfig,'PaperPositionMode','auto');
            print(netfig,'-depsc',[pathstr slash 'Figures' slash name '.eps']);
            
            print(netfig,'-dtiff',[pathstr slash 'Figures' slash name '.tif']);
            close(netfig);
            drawnow;
             multiWaitbar('Saving results and making summary figure',1,'Color','g');
            end
        end
        data = processed_analysis(j);
         [pathstr, name] = fileparts(processed_analysis(j).filename);
        save([pathstr slash 'analysis-' name '.mat'],'data');
        
       
        
        pause(1);
        multiWaitbar('CloseAll','Name','');
    end
    
    % Finally, save it in the same location as analysis.mat
    if(~ispc)
        TARGET_PATH = strrep(TARGET_PATH,'\','/');
    end
    cprintf('*blue','%s\n', ['Saving to ' TARGET_PATH Folder slash 'processed_analysis.mat']);
    save([TARGET_PATH Folder slash 'processed_analysis.mat'],'processed_analysis');
    
end
% matlabpool close