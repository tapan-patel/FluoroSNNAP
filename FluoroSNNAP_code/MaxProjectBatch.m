function MaxProjectBatch(varargin)
if(~matlabpool('size'))
    matlabpool open
end
if(nargin>0)
    FNames = varargin{1};
    SOURCE_PATH = [];
else
SOURCE_PATH = Read_Params('SOURCE_PATH');
FNames = Read_FNames;
end
for i=1:length(FNames)
    Folder = FNames{i};
  
    files = dir([SOURCE_PATH Folder '/*.tif']);
    
     % Exclude files that begin with MAX or Segmentation or end with -file
     exclude_idx = [];
     for j=1:numel(files)
         if(~isempty(strfind(files(j).name,'Segmentation')) || ...
             ~isempty(strfind(files(j).name, 'MAX')) || ...
             ~isempty(strfind(files(j).name, '-file')))
             exclude_idx = [exclude_idx j];
         end
     end
     files(exclude_idx) = [];
     
     % Loop through each file and compute max projection image - to be used
     % with segmentation
     parfor j=1:numel(files)
         MaxProjection([SOURCE_PATH Folder '/'],files(j).name);
     end
end