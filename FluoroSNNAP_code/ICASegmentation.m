function Lica = ICASegmentation(fn,flims,nPCs,PCuse,mu,nIC,smwidth,thresh,arealim,plotting,dsamp,badframes,outputdir)
multiWaitbar('CloseAll','Name','');
[~,y] = fileparts(fn);
multiWaitbar('Performing PCA',0,'Name',y,'CancelFcn', @(a,b) disp( ['Cancel ',a] ) );
multiWaitbar('Performing ICA',0,'Color','b','CancelFcn', @(a,b) disp( ['Cancel ',a] ) );
multiWaitbar('Performing segmentation',0,'Color','r','CancelFcn', @(a,b) disp( ['Cancel ',a] ) );

[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs,dsamp,outputdir,badframes);
multiWaitbar('Performing PCA',1);

[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals,PCuse,mu,nIC);
multiWaitbar('Performing ICA',1,'Color','b');

multiWaitbar('Performing segmentation',0,'Color','r','Busy');
[ica_segments, segmentlabel, segcentroid,~,Lbig,ica_filtersbw] = CellsortSegmentation(ica_filters, smwidth, thresh, arealim, plotting);
multiWaitbar('Performing segmentation',.8,'Color','r');
nsegments = size(ica_segments,1);
Lica_stack = cell(nsegments,1);
for i=1:nsegments
    Lica_stack{i} = logical(squeeze(ica_segments(i,:,:)));
end
% Remove ROIs where one ROI is completely within another ROI.
% Do not use bwlabel because ROIs that are touching each other with be
% labeled as a single ROI. Need to manually renumber ROIs.

remove = [];
% Find indices of all ROIs
IDX = cell(nsegments,1);
for i=1:nsegments
    IDX{i} = find(Lica_stack{i});
end
for i=1:nsegments
    for j=i+1:nsegments
        idx1 = IDX{i};
        idx2 = IDX{j};
        n1 = length(idx1); n2 = length(idx2);
        if(n1==0)
            remove = [remove i];
        
        elseif(n2==0)
            remove = [remove j];
        end
         if(n1>n2 && length(intersect(idx1,idx2))/n2>0.9)
            remove = [remove j];
        elseif(n2>n1 && length(intersect(idx2,idx1))/n1>0.9)
            remove = [remove i];
        end
    end
end
remove = unique(remove);
Lica_stack(remove) = [];

nsegments = size(Lica_stack,1);
I=imread(fn);
Lica = zeros(size(I));

for i=1:nsegments
    Lica(Lica_stack{i}) = i;
end
% make sure the ROIs fall between area limits
C = regionprops(Lica,'Area');
A = [C.Area];
for i=1:max(Lica(:))
    if(A(i)<arealim(1) || A(i)>arealim(2))
        Lica(Lica==i)=0;
    end
end
% Renumber ROIs
L = Lica;
for k=1:3
    nROIs = max(L(:));
    for j=1:nROIs
        if(isempty(find(L==j,1))) % This ROI does not exist
            for i=j:nROIs-1
                L(L==i+1) = i;
            end
            nROIs = nROIs-1;
        end
    end
end
[folder,file] = fileparts(fn);
ica = 1;
savename = [folder '/Segmentation-' file '.mat'];
save(savename,'L','ica');
disp(['Segmentation saved to ' savename]);
multiWaitbar('Performing segmentation',1,'Color','r');
multiWaitbar('CloseAll');