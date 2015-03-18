function MaxProjection(path,imgname)
% Read in a tiff stack and generate max projection image
disp(['Processing ' path imgname]);
%% Get file name
info = imfinfo([path imgname]);
I = zeros(info(1).Height,info(1).Width);
N = numel(info);
%%
for i=1:numel(info)
    if(~mod(i,500))
        disp([num2str(i) ' frames of ' num2str(N) ' read.']);
    end
    img = double(imread([path imgname],'Info',info,'Index',i));
    max_idx = img>I;
    I(max_idx) = img(max_idx);
end

%% Save 
% figure;
% imshow(imadjust(uint8(I/numel(info))),[]);
% title([path imgname]);
imwrite(uint16(I),[path 'MAX-' imgname],'tif');
disp(['Max projection image writen to ' path 'MAX-' imgname]);