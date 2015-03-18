%%
fn = 'baseline.tif';
flims = [1 1200];
nPCs = 150;
dsamp = [];
badframes = [];
outputdir = [];
[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(fn, flims, nPCs,dsamp,outputdir,badframes);
PCuse = 1:100;
mu = 1;
nIC = 100;
smwidth = 4;
thresh = 2;
arealim = [50 500];
plotting = 1;
%%
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals,PCuse,mu,nIC);

%%
[ica_segments, segmentlabel, segcentroid,~,Lbig,ica_filtersbw] = CellsortSegmentation(ica_filters, smwidth, thresh, arealim, plotting);
%%
cell_sig = CellsortApplyFilter(fn, ica_segments);
%%
Lnew = zeros(size(Imean));
for i=1:122
    bw = squeeze(ica_segments(i,:,:))~=0;
    Lnew(bw) = i;
end
    
%%
figure; imshow(imadjust(uint16(Imean)));
hold on
for k=1:116
B = bwboundaries(L==k);

for i=1:length(B)
    plot(B{i}(:,2),B{i}(:,1),'r');
end
end
%%
L = zeros(size(Imean));
for i=1:100
    L = L + squeeze(ica_filtersbw(:,:,i));
end
%%
figure; imshow(imadjust(uint16(Imean)));
hold on
for i=1:size(segcentroid,1)
    plot(segcentroid(i,1),segcentroid(i,2),'ro');
end