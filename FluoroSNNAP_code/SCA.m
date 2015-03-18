function Cluster = SCA(s,varargin)
% Synchronization cluster analysis: Determines the number clusters and the
% participation index.
% Input: MxT matrix of time-series fluorescence data
% Ref: Synchronization measurement of multiple neuronal populations, Li
% 2007 J Neurophysiol 98: 3341-3348
% Ref: sclae-free topology of the CA2 hippocampal network, Li 2010
% Biophysical Journal vol 99, 1733-1741

save_flag = 0;
DISPLAY = 0;
F = s.F_cell;
wb = 0; % show wait bar
analysis_type = 'phase'; % Default to event synchronization from Li 2010,
% other options are circular variance of phase or equal-time correlation
lag = .1*s.fps; % events occuring closer than lag are considered simulatenous
try
    load('params.mat');
catch
    errordlg('Could not load params.mat. Please run Analysis -> Preferences -> Revert to default','File not found','modal');
end

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

if(wb)
    msg = sprintf('Performing synchronization cluster analysis.\nThis may take several minutes.\nThis message will self-destruct');
    wbh = msgbox(msg,'Please wait');
    delete(findobj(wbh,'string','OK')); drawnow;
end
% Preprocessing
[M,T] = size(F);
% if(~isfield(s,'SynchroCluster') || isempty(s.SynchroCluster))
% Need to compute the firingr rate for surrogate testing
r = zeros(M,T); % Firing rate
C = zeros(M,M); % Event synchronization matrix
if(params.sca_type==1)
    analysis_type = 'phase';
elseif(params.sca_type==2)
    analysis_type = 'correlation';
elseif(params.sca_type==3)
    analysis_type = 'entropy';
end
switch analysis_type
    case 'event'
        
        for i=1:M
            for j=i:M
                c_ij = SimultaneousEvents(s.Spikes_cell{i},s.Spikes_cell{j},lag);
                c_ji = SimultaneousEvents(s.Spikes_cell{j},s.Spikes_cell{i},lag);
                mi = length(s.Spikes_cell{i});
                mj = length(s.Spikes_cell{j});
                if(mi*mj>0)
                    C(i,j) = (c_ij+c_ji)/sqrt(mi*mj);
                    C(j,i) = C(i,j);
                end
            end
        end
        
        t = 0:1/s.fps:T/s.fps-1/s.fps;
        for i=1:M
            fr = FiringRate(s.Spikes_cell{i},s.fps,T);
            r(i,:) = fr;
        end
    case 'correlation'
        
        Z = zscore(F'); Z = Z';
        
        %% Equal-time correlation matrix & Eigenvalue decomp
        
        
        for i=1:M
            for j=i:M
                C(i,j) = 1/T*Z(i,:)*Z(j,:)'; % Simple dot product
                C(j,i) = C(i,j); % Symmetric matrix with diag =1
            end
        end
    case 'phase'
        % Recompute phase everytime this function is called since user may
        % have updated spike with the CaGUI so just to be sure, recalculate
        % instantaneous phase
        s.phase = [];
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
        phi = s.phase;
        Cluster.phase = phi;
        for i=1:M
            for j=i:M
                if(isempty(s.Spikes_cell{i}) || isempty(s.Spikes_cell{j}))
                    C(i,j) = 0;
                    C(j,i) = 0;
                else
                    deltaphi = mod(phi(i,:)-phi(j,:),2*pi);
                    C(i,j) = sqrt(mean(cos(deltaphi))^2+mean(sin(deltaphi))^2);
                    C(j,i) = C(i,j);
                end
            end
        end
    case 'entropy'
        if(~isfield(s,'phase'))
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
        phi = s.phase;
        for i=1:M
            for j=i:M
                if(isempty(s.Spikes_cell{i}) || isempty(s.Spikes_cell{j}))
                    C(i,j) = 0;
                    C(j,i) = 0;
                else
                    [c,lags] = xcorr(mod(phi(i,:),2*pi),mod(phi(j,:),2*pi));
                    [~,idx] = max(c);
                    tau = lags(idx);
                    a = mod(phi(i,:),2*pi);
                    b = mod(phi(j,:),2*pi);
                    
                    if(tau<0)
                        if(tau<-200)
                            tau = -1;
                        end
                        b = b(-1*tau:end);
                        a = a(1:length(b));
                        
                    elseif(tau>0)
                        if(tau>200)
                            tau = 1;
                        end
                        a = a(tau:end);
                        b = b(1:length(a));
                    end
                    [I,HA,HB] = MutualInformation(a,b);
                    tmp = I/(max([HA HB]));
                    if(isnan(tmp) || isinf(tmp))
                        C(i,j) = 0;
                        C(j,i) = 0;
                    else
                        C(i,j) = tmp;
                        C(j,i) = C(i,j);
                    end
                end
            end
        end
end

% C is a real symmetric matrix, with real positive eigenvalues. If the time
% series are fully uncorrelated, all eigenvalues will distribute around 1.
% If all the series are fully correlated, max eigenval = M and the rest
% fall to 0.
[V,D] = eig(C);
D = diag(D);
%% Significance testing with surrogates. Construct surrogate time series
%% using amplitude-adjusted fourier transform (MATs toolbox) and perform
%% equal-time correlation analysis N times to give a mean and std for
%% surrogate eigenvalues. For event-synchronization, use the actual
%% spike train to compute time dependent firing rate and use that to
%% generate surrogate spikes from inhomogeneous Poisson process
N = params.sca_N; % Create surrogates 20 times
Spikes_cell = s.Spikes_cell;

Dsur = zeros(M,N); % Columwise eigenvalues for surrogate channels

% Generate surrogate spikes for event-synchronization using
% inhomogenous poisson process for total N times - this is the longest
% part of the code

SurSpikes = cell(M,N);
% Generate a Poisson spike train
if(strcmp(analysis_type,'event'))
    for n=1:N
        if(params.parallel)
            parfor i=1:M
                x = InHomoPoisSpkGen(r(i,:));
                SurSpikes(i,n) = {x};
            end
        else
            for i=1:M
                x = InHomoPoisSpkGen(r(i,:));
                SurSpikes(i,n) = {x};
            end
        end
    end
end
Rall = zeros(M,M,N);
for n=1:N % Index for surrogate testing
    multiWaitbar('Synchronization Cluster Analysis',n/N);
    switch analysis_type
        case 'event'
            R = zeros(M,M);
            % Compute surrogate correlation matrix R
            for i=1:M
                for j=i:M
                    r_ij = SimultaneousEvents(SurSpikes{i,n},SurSpikes{j,n},lag);
                    r_ji = SimultaneousEvents(SurSpikes{j,n},SurSpikes{i,n},lag);
                    mi = length(SurSpikes{i});
                    mj = length(SurSpikes{j});
                    if(mi*mj>0)
                        R(i,j) = (r_ij+r_ji)/sqrt(mi*mj);
                        R(j,i) = R(i,j);
                    end
                end
            end
            % Make sure this is a symmetric matrix
            if(~isequal(R,R'))
                error('Surrogate matrix is not symmetric');
            end
        case 'correlation'
            Fsur = zeros(M,length(F(1,:)));
            for i=1:M
                Fsur(i,:) = AAFTsur(F(i,:),1); % AAFT surrogate
                %           Fsur(i,:) = F(i,randperm(T+1));  % Random permutation surrogate
            end
            Fsur = zscore(Fsur'); Fsur = Fsur';
            R = zeros(M,M);
            for i=1:M
                for j=i:M
                    R(i,j) = 1/T*Fsur(i,:)*Fsur(j,:)';
                    R(j,i) = R(i,j);
                end
            end
        case 'phase'
            phi = s.phase;
            for i=1:M
                for j=i:M
                    if(isempty(s.Spikes_cell{i}) || isempty(s.Spikes_cell{j}))
                        R(i,j) = 0;
                        R(j,i) = 0;
                    else
                        phi_j = AAFTsur(phi(j,:));
                        deltaphi = mod(phi(i,:)-phi_j',2*pi);
                        R(i,j) = sqrt(mean(cos(deltaphi))^2+mean(sin(deltaphi))^2);
                        R(j,i) = R(i,j);
                    end
                end
            end
        case 'entropy'
            phi = s.phase;
            for i=1:M
                for j=i:M
                    phi_j = AAFTsur(phi(j,:));
                    [c,lags] = xcorr(mod(phi(i,:),2*pi),mod(phi_j,2*pi));
                    [~,idx] = max(c);
                    tau = lags(idx);
                    a = mod(phi(i,:),2*pi);
                    b = mod(phi_j,2*pi);
                    if(tau<0)
                        b = b(-1*tau:end);
                        a = a(1:length(b));
                        
                    elseif(tau>0)
                        a = a(tau:end);
                        b = b(1:length(a));
                    end
                    [I,HA,HB] = MutualInformation(a,b);
                    tmp = I/(max([HA HB]));
                    if(isnan(tmp) || isinf(tmp))
                        R(i,j) = 0;
                        R(j,i) = 0;
                    else
                        R(i,j) = tmp;
                        R(j,i) = R(i,j);
                    end
                end
            end
    end
    Rall(:,:,n) = R;
    [~,Dtmp] = eig(R);
    Dsur(:,n) = diag(Dtmp);
end

lambda_sur = mean(Dsur,2); % Mean surrogate eigenvalues
SD_sur = std(Dsur,0,2);

SI = zeros(M,1);

for i=1:M
    if(D(i)>(lambda_sur(i)+2*SD_sur(i)))
        SI(i) = (D(i)-lambda_sur(i))/(M-lambda_sur(i));
    end
end

%% Compute the participation index for each cluster. Ignore SI<0.01

Nclusters = sum(SI>=0.01);
PI = zeros(M,Nclusters);
cluster_i = find(SI>=0.01);

for k=1:Nclusters
    lambda = D(cluster_i(k));
    vk = V(:,cluster_i(k));
    
    PI(:,k) =lambda*vk.^2;
end
% There is a minimum number of neurons required to create a synchronization
% cluster. User specified this value
[~,idx] = max(PI,[],2);
Cluster_size = zeros(size(PI,2),1);
for i=1:length(Cluster_size)
    Cluster_size(i) = nnz(idx==i);
end
Nclusters = nnz(Cluster_size>params.sca_size);
idx = find(Cluster_size<=params.sca_size);
PI(:,idx) = [];
% Neurons whose PI>0.05 are part of an assembly. Single neuron can be
% part of >1 assembly -> hub neurons.
assembly_contents = cell(size(Nclusters));
for k=1:Nclusters
    assembly_contents(k) = {find(PI(:,k)>=0.05)};
end

Cluster.SI = SI;
Cluster.Num = Nclusters;
Cluster.PI = PI;
Cluster.EigenVal = D;
Cluster.Surrogate = lambda_sur;
Cluster.C = C;

Cluster.Assembly_Contents = assembly_contents;

if(wb)
    close(wbh)
end
%% Make figure and save if requested
if(DISPLAY || save_flag)
    save_directory = pwd;
    if(DISPLAY)
        h = figure;
    else
        h = figure('visible','off');
    end
    % Plot the size of each cluster
    ClusterSize = zeros(length(Cluster.Assembly_Contents),2);
    
    for i=1:Cluster.Num
        ClusterSize(i,1) = i;
        ClusterSize(i,2) = length(Cluster.Assembly_Contents{i});
    end
    stem(ClusterSize(:,1),ClusterSize(:,2));
    xlim([0 Cluster.Num+1]);
    xlabel('Cluster Number'); ylabel('# of cells belonging to a cluster');
    legend(num2str(Cluster.SI(Cluster.SI>.01)))
    title(['SCA for ' save_directory]);
    if(save_flag)
        print(h,'-depsc',save_directory(1:end-4));
    end
end

function c=SimultaneousEvents(x,y,lag)
% Given two spike trains x and y, the number of times that an event appears
% in x shortly after it appears in y within x specified time lag
c = 0;

mx = length(x); % Number of spikes in x
my = length(y);  % Number of spikes in y

% For x given spike at t=t0 in y, figure out if x spike occured in x within
% t=t0 to t=t0+lag. If the spike occured at t=t0 in both x and y, c=c+.5,
% else c= c+1
for i=1:my
    for j=1:mx
        if(x(j)-y(i) <=lag &  x(j)-y(i) >0)
            c = c+1;
        elseif(x(j)==y(i))
            c = c+.5;
        end
    end
end