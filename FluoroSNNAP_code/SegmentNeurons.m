function SegmentNeurons(varargin)

global filename BW bkg_coords;
if(nargin>0)
    filename = varargin{1};
else
    [imgname, path] = uigetfile('*.tif','Select an image to segment');
    filename = [path '/' imgname];
end
I = imread(filename);
roiwindow = CROIEditor(I); drawnow;

while(roiwindow.isvalid & isempty(roiwindow.number))
    pause(1)
    
    if(roiwindow.isvalid & roiwindow.number)
        break
    end
end
if(roiwindow.isvalid)
    [BW,bkg_coords] = roiwindow.getROIData;
    [folder, imgname] = fileparts(filename);
    % If max projection was done, need to remove "MAX-" from the file.
    imgname = strtok(imgname,'MAX-');
    save([folder '/Segmentation-' imgname],'BW','bkg_coords');
    imwrite(BW,[folder '/Segmentation-' imgname '.tif']);
    disp(['Segmentation saved to ' folder '/Segmentation-' imgname])
    delete(roiwindow);
end