function [A,rise_time,fall_time,CV] = GetTransients(s,celltype)
warning('off')
% Extract snippets of transients and compute the amplitudes
% If neurons and astrocytes have been discriminated using CellTypeNMDA
% function, use that information to compute transients from neurons only.
% Else use all cells

[~,fname] = fileparts(s.filename);
N = size(s.F_cell,1);
if(nargin==1)
    neurons = true(N,1);
else
    neurons = celltype==2;
end

s.F_cell = s.F_cell(neurons,:);
s.Spikes_cell = s.Spikes_cell(neurons);

[N,~] = size(s.F_cell);
Transient = cell(N,1);
A = zeros(N,1);
CV = zeros(N,1);
rise_time = cell(N,1);
fall_time = cell(N,1);
for i=1:N
    multiWaitbar('Computing calcium event transient kinetics',i/N,'Name',fname);
    F = s.dF_cell(i,:);
    Y = F;
    %     Y = wden(F,'heursure','s','one',5,'sym8');
    spikes = s.Spikes_cell{i};
    % spikes is a vector of timestamps of the onset of transient.
    % Store the transient begining at he spike onset and
    if(isempty(spikes))
        A(i) = 0; rise_time{i} = 0; fall_time{i} = 0; CV(i) = 0;
    else
        amplitudes = nan(1,length(spikes));
        rise = nan(1,length(spikes));
        fall = nan(1,length(spikes));
        %     transients = zeros(length(spikes),101);
        spikes(end+1) = size(s.F_cell,2); % Add artificial spike at the end but loop to spikes-1
        
        for j=1:length(spikes)-1
            try
                G = Y(spikes(j)-3:min([spikes(j)+5*s.fps spikes(j+1)-3]));
                
                amplitudes(j) = range(G);
            end
            % Get rise and fall time for each transient
            try
                
                [rise(j), fall(j)] = TransientKinetics(G,s.fps);
            end
        end
        
        %     transients(amplitudes==0,:) = [];
        try
            A(i) = mean(amplitudes(~isnan(amplitudes)));
            CV(i) = std(amplitudes(~isnan(amplitudes)))/A(i);
            %     Transient(i) = {transients};
        end
        rise_time{i} = rise(~isnan(rise));
        fall_time{i} = fall(~isnan(fall));
    end
end

% Plot the snippets
% figure
% t = 0:1/s.fps:100/s.fps;
% for i=1:N
%     tr = Transient{i};
%     for j=1:size(tr,1);
%         plot(t,tr(j,:));
%         hold on
%     end
% end

% Plot the distribution of amplitudes
% figure
% for i=1:N
%     try
%         plot(i,A{i},'b.');
%         plot(i,mean(A{i}),'r.')
%         CV(i) = std(A{i})/mean(A{i});
%         hold on
%     end
% end
% xlabel('Neuron ID');
% ylabel('Amplitude per transient');
% ci = bootci(100,@(x) mean(x),CV);
% title(['95% CI for CV = [' num2str(ci(1)) ' ' num2str(ci(2)) ']']);
% compiled_amplitudes = [];
% for i=1:N
%     compiled_amplitudes = [compiled_amplitudes A{i}];
% end
%
% figure
% hist(compiled_amplitudes);
% xlim([0 .8]);
% xlabel('Mean amplitude per neuron'); ylabel('Frequency');
% title(['Mean = ' num2str(mean(compiled_amplitudes)) '. Std = ' num2str(std(compiled_amplitudes)/sqrt(numel(compiled_amplitudes)))]);