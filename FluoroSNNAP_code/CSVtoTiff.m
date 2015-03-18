function CSVtoTiff(filename)
% Input: filename of a .csv file
% Input file is a comma-separated file; columns are the ROIs
% (neurons), rows are fluorescence intensities at each frame. The number of
% rows = number of frames. The number of columns = number of ROIs

% Since the CalciumAnalysis software operates on tiff image stack, we need
% to convert the CSV file to "mock" image stack.
M = dlmread(filename);
[T,N] = size(M); % T=number of frames, N = number of neurons (ROIs)

[folder,file,ext] = fileparts(filename);
r = max([600 2*N]);
Istack = zeros(r,r,N,'uint16');

% Place the N neurons on a 400x400 image. Generate centroid positions.
x = randsample(r,N);
y = randsample(r,N);

% Create a label matrix with these centroids
L = zeros(r);
for i=1:N
    L(x(i),y(i)) = i;
end
L=imdilate(L,strel('disk',5));
%Save the segmentation
save([fullfile(folder,['Segmentation-' file]) '.mat'],'L');
hwait = msgbox(sprintf('Please wait. Converting .csv to .tif\n%s',filename));
% Now loop through the frames and place the frame intensity from the csv
% file of each ROI in the image stack
for i=1:T
    for j=1:N
        Istack(x(j),y(j),i) = uint16(M(i,j));
    end
    Istack(:,:,i) = imdilate(Istack(:,:,i),strel('disk',5));
end

imwrite(Istack(:,:,1),[fullfile(folder,file) '.tif']);
for i=2:T
    imwrite(Istack(:,:,i),[fullfile(folder,file) '.tif'],'WriteMode','append');
end
try
    delete(hwait);
end
