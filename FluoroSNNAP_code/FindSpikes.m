function TE = FindSpikes(F,varargin)
% %% Detect spikes using a moving 5 point line fit
% warning off
% % Y = zscore(F); % Observed fluorescence
% Y = wden(F,'heursure','s','one',5,'sym8');
% % Y = smooth(Y,'sgolay'); Y = Y';
% N = numel(Y);
% n = 10;
% fps=20;
% x = (0:n-1)/fps;
thr = -.1;
% Ns = 20;
if nargin > 1
    for i = 1:2:length(varargin)-1
        if isnumeric(varargin{i+1})
            eval([varargin{i} '= [' num2str(varargin{i+1}) '];']);
        else
            eval([varargin{i} '=' char(39) varargin{i+1} char(39) ';']);
        end
    end
end
Signal = F;
SFr = 20;
Wid = [0.5 1.0];
Ns = 20; 
option = 'l'; L=thr; wname = 'sym2'; PltFlg = 1; CmtFlg =1;

wfam = {'bior1.5','bior1.3','sym2','db2','haar'};

if sum(strcmp(wname,wfam)) == 0
    error('unknown wavelet family')
elseif CmtFlg == 1
    disp(['wavelet family: ' wname])
    to = clock;
end

%make sure signal is zero-mean
Signal = Signal - mean(Signal);

Nt = length(Signal);      %# of time samples

%define relevant scales for detection
W = determine_scales(wname,Wid,SFr,Ns);

%initialize the matrix of thresholded coefficients
ct = zeros(Ns,Nt);

%get all coefficients 
c = cwt(Signal,W,wname);  

%define detection parameter
Lmax = 36.7368;       %log(Lcom/Lom), where the ratio is the maximum 
                    %allowed by the current machine precision
L = L * Lmax;

%initialize the vector of spike indicators, 0-no spike, 1-spike
Io = zeros(1,Nt);

%loop over scales
for i = 1:Ns
    
    %take only coefficients that are independent (W(i) apart) for median
    %standard deviation
    
    Sigmaj = median(abs(c(i,1:round(W(i)):end) - mean(c(i,:))))/0.6745;
    Thj = Sigmaj * sqrt(2 * log(Nt));     %hard threshold
    index = find(abs(c(i,:)) > Thj);
    if isempty(index) & strcmp(num2str(option),'c')
        %do nothing ct=[0];
    elseif isempty(index) & strcmp(num2str(option),'l')
        Mj = Thj;
        %assume at least one spike
        PS = 1/Nt;
        PN = 1 - PS;
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];    %decision threshold
        DTh = abs(DTh) * (DTh >= 0);                 %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        if isempty(ind)
            %do nothing ct=[0];
        else
            ct(i,ind) = c(i,ind);
        end
    else
        Mj = mean(abs(c(i,index)));       %mean of the signal coefficients
        PS = length(index)/Nt;            %prior of spikes
        PN = 1 - PS;                        %prior of noise
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];   %decision threshold
        DTh = abs(DTh) * (DTh >= 0);         %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        ct(i,ind) = c(i,ind);
    end
    
    %find which coefficients are non-zero
    Index = ct(i,:) ~= 0;
    
    %make a union with coefficients from previous scales
    Index = or(Io,Index);
    Io = Index;
end

TE = parse(Index,SFr,Wid);

if PltFlg == 1
    close all
    figure(1)
    scale = 64./[max(abs(c),[],2) * ones(1,Nt)];
    temp = zeros(1,Nt);
    temp(TE) = 1;
    image(flipud(abs(c)) .* scale)
    colormap pink
    ylabel('Scales')
    Wt = [fliplr(W)];
    set(gca,'YTick',1:Ns,'YTickLabel',Wt,'Position',[0.1 0.2 0.8 0.6], ...
        'XTick',[])
    title(['|C| across scales: ' num2str(W)])
    ah2 = axes;
    set(ah2,'Position',[0.1 0.1 0.8 0.1])
    plot(temp,'o-m','MarkerSize',4,'MarkerFaceColor','m')
    set(gca,'YTick',[],'XLim',[1 Nt])
    xlabel('Time (samples)')
    ylabel('Spikes')
    
    figure(2)
    plot(Signal,'Color',[0.7 0.7 0.7],'LineWidth',2)
    hold on
    plot(ct','-o','LineWidth',1,'MarkerFaceColor','k', ...
        'MarkerSize',4)
    xlabel('Time (samples)')
    ylabel('Coefficients')
    set(gca,'XLim',[1 Nt])
end

if CmtFlg == 1
    disp([num2str(length(TE)) ' spikes found'])
    disp(['elapsed time: ' num2str(etime(clock,to))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Scale = determine_scales(wname,Wid,SFr,Ns)
  
%Ns - # of scales  

dt = 1/SFr;  %[msec]

%signal sampled @ 1 KHz  
Signal = zeros(1,1000);
%create Dirac function
Signal(500) = 1;
  
Width = linspace(Wid(1),Wid(2),Ns);

%infinitesimally small number
Eps = 10^(-15);

ScaleMax = 3;
ScaleMax = ScaleMax*SFr;

switch num2str(wname)
  
 case 'haar'
  for i = 1:Ns
    Scale(i) = Width(i)/dt - 1; 
  end
 case 'db2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'sym2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of positive slope zero crossings
    IndZeroCross = find(IndDer == 1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.3'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
   for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.5'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 otherwise
  error('unknown wavelet family')
end

NaNInd = isnan(Scale);

if sum(NaNInd) > 0
  warning(['Your choice of Wid is not valid given' ...
        ' the sampling rate and wavelet family'])
  if NaNInd(1) == 1
    disp(['Most likely Wid(1) is too small'])
  elseif NaNInd(Ns) == 1
    disp(['Most likely Wid(2) is too large'])
    disp(['Change the value on line: ''ScaleMax = 2'' to something larger'])
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn = parse(Index,SFr,Wid);

%This is a special function, it takes the vector Index which has 
%the structure [0 0 0 1 1 1 0 ... 0 1 0 ... 0]. This vector was obtained
%by coincidence detection of certain events (lower and upper threshold
%crossing for threshold detection, and the appearance of coefficients at
%different scales for wavelet detection). 
%The real challenge here is to merge multiple 1's that belong to the same
%spike into one event and to locate that event

Refract = 1.5 * Wid(2);    %[ms] the refractory period -- can't resolve spikes 
                           %that are closer than Refract;
Refract = round(Refract * SFr);

Merge = mean(Wid);      %[ms] merge spikes that are closer than Merge, since 
                        %it is likely they belong to the same spike

Merge = round(Merge * SFr);   


Index([1 end]) = 0;   %discard spikes located at the first and last samples

ind_ones = find(Index == 1);    %find where the ones are

if isempty(ind_ones)
    TE = [];
else
    temp = diff(Index);  %there will be 1 followed by -1 for each spike
    N_sp = sum(temp == 1); %nominal number of spikes
    
    lead_t = find(temp == 1);  %index of the beginning of a spike
    lag_t = find(temp == -1);  %index of the end of the spike
    
    for i = 1:N_sp
        tE(i) = ceil(mean([lead_t(i) lag_t(i)]));
    end
   
    i = 1;        %initialize counter
    while 0 < 1
        if i > (length(tE) - 1)
            break;
        else
            Diff = tE(i+1) - tE(i);
            if Diff < Refract & Diff > Merge
                tE(i+1) = [];      %discard spike too close to its predecessor
            elseif Diff <= Merge
                tE(i) = ceil(mean([tE(i) tE(i+1)]));  %merge
                tE(i+1) = [];                         %discard
            else
                i = i+1;
            end
        end
    end 
    TE = tE;
end

fcn = TE;

% 
% slope = @(x,y) ( x*y' - 1/n*sum(x)*sum(y))/( sum(x.^2) - 1/n*sum(x)^2);
% 
% breakpts = 1:n/2:N-n;
% M = zeros(numel(breakpts),1);
% 
% for i=1:numel(breakpts)
%     y = Y(breakpts(i):breakpts(i)+n-1);
%     M(i) = slope(x,y);
% end
% 
% % Threshold the slopes - this is tricky. Threshold = average of absolute
% % slopes + n*std(M). Set n=2. This should work for most noisy signals. But,
% % if there is a huge spike in the signal, figure out how many sd above mean
% % this corresponds to (=k) and set n= avg(2,k)
% 
% % L=200;
% % d=100;
% % thresh = zeros(numel(1:d:length(M)-L-d)+1,1);
% % cntr = 1;
% % try
% %     for i=1:d:length(M)-L-d
% %
% %         nsd = max(M(i:i+L))/std(M(i:i+L));
% %         thresh(cntr) = nsd;
% %         cntr = cntr+1;
% %     end
% %     % Take care of the omitted samples
% %     if(length(M)-L-d>0)
% %         nsd = max(M(length(M)-L-d:end))/std(M(length(M)-L-d:end));
% %         thresh(end) = nsd;
% %     else
% %         thresh(end) = thresh(end-1);
% %     end
% % end
% % threshold = mean(thresh)+std(M);
% % threshold = threshold/min(Y);
% % threshold = max([threshold .06]);
% % indx = 1:numel(breakpts);
% %
% % spikes_ind = indx(M/min(Y)>threshold);
% 
% 
% % Detect spikes by Continuous wavelet transform from scales 2-20 and 20-60.
% % 2-20 will detect somewhat short lived transients (normal). 20-60 scales
% % will pick up long-lived transients, typical of injured cells
% 
% wt = cwt(M,1:Ns,'db4');
% % [c,l] = wavedec(M',5,'db4');
% % sigma = wnoisest(c,l,1);alpha=2;
% 
% % thr = min([wbmpen(c,l,sigma,alpha) 2]);
% 
% 
% % thr = max([max(wt(:))/2 1]);
% % thr = max([mean(max(wt,[],2))/2 .5]);
% 
% % Loop through the scales and find the local maxima above 20
% spikes_ind = zeros(1,length(M));
% for i=1:Ns
%     [~,val] = findpeaks(wt(i,:),'minpeakheight',0.02);
%     
%     spikes_ind(val(Y(breakpts(val))>thr)) = 1;
% % spikes_ind(val(wt(i,val)>thr)) = 1;
%   
%     %     spikes_ind(val) = 1;
%     %     [val,locs] = findpeaks(wt(:,i));
%     %     if(~isempty(locs))
%     %         if(max(val)>=10 & M(i)>=5)
%     %             spikes_ind(i) = i;
%     %         end
%     %     end
% end
% spikes_ind = find(spikes_ind);
% if(~isempty(spikes_ind))
%     spikes_ind = sort(spikes_ind,'descend');
%     spikes_ind([3 diff(spikes_ind)]~=-1);
% end
% 
% % Now look at larger scales
% spikes_ind2 = zeros(1,length(M));
% % wt = cwt(M,20:80,'db4');
% % thr = max([mean(max(wt,[],2)) .5]);
% % for i=1:size(wt,1)
% %     [~,val] = findpeaks(wt(i,:));
% % 
% %     spikes_ind2(val(wt(i,val)>thr & M(val)'>.1)) = 1;
% %   
% %     %     spikes_ind(val) = 1;
% %     %     [val,locs] = findpeaks(wt(:,i));
% %     %     if(~isempty(locs))
% %     %         if(max(val)>=10 & M(i)>=5)
% %     %             spikes_ind(i) = i;
% %     %         end
% %     %     end
% % end
% % 
% % spikes_ind2 = find(spikes_ind2);
% % spikes_ind2 = sort(spikes_ind2,'descend');
% % % Remove consecutive spikes
% % if(~isempty(spikes_ind2))
% %     spikes_ind2 = spikes_ind2([3 diff(spikes_ind2)]~=-1);
% % end
% % 
% % 
% % 
% % spikes_ind3 = zeros(1,length(M));
% % wt = cwt(M,20:80,'db4');
% % thr = max([mean(max(wt(:,floor(length(M/2))),[],2)) .5]);
% % for i=1:size(wt,1)
% %     [~,val] = findpeaks(wt(i,1:floor(length(M/2))));
% % 
% %     spikes_ind3(val(wt(i,val)>thr & M(val)'>.1)) = 1;
% %   
% %     %     spikes_ind(val) = 1;
% %     %     [val,locs] = findpeaks(wt(:,i));
% %     %     if(~isempty(locs))
% %     %         if(max(val)>=10 & M(i)>=5)
% %     %             spikes_ind(i) = i;
% %     %         end
% %     %     end
% % end
% % 
% % spikes_ind3 = find(spikes_ind3);
% % spikes_ind3 = sort(spikes_ind3,'descend');
% % % Remove consecutive spikes
% % if(~isempty(spikes_ind3))
% %     spikes_ind3 = spikes_ind3([3 diff(spikes_ind3)]~=-1);
% % end
% %     % If there is a gap of 2 or less, merge them
% spikes_total = [spikes_ind];
% spikes_total = (sort(unique(spikes_total),'descend'));
% 
% if(~isempty(spikes_total))
%     spikes_total = spikes_total([3 diff(spikes_total)]~=-1);
%     spikes = sort(breakpts(spikes_total));
% else
%     spikes = [];
% end
% % spikes = spikes_total;
