function Data = AnalyzeData(varargin)

% Reads list of folders that contain .tif image stacks to compute the
% time-series fluorescence. All .tif files in the given folder will be
% processed, unless indicated under the exclusion list (Read_FNames.m).
% Options are to compute full-frame fluorescence or individual cell
% fluorescence or both. If individual cell fluorescence is requested, a
% segmentation file must be present in the folder that contains the
% original image stack.

% Convention ** Segmentation file for a stack named "filename.tif" must be
% named "Segmentation-filename.tif"

% Output is saved as analysis.mat file in the same folder. Run PostProcess
% function to detect spikes and other analyses.


if(nargin>0)
    SOURCE_PATH = [];
    TARGET_PATH = [];
    ANALYSIS_TYPE = 'BOTH';
    SAVE_RESULTS = 1;
    FNames = varargin{1};
    if(nargin==2)
        files = varargin{2};
        files.name = files;
    end
else
    %% Get the params
    SOURCE_PATH = Read_Params('SOURCE_PATH');
    TARGET_PATH = Read_Params('TARGET_PATH');
    ANALYSIS_TYPE = Read_Params('ANALYSIS_TYPE'); % Full frame, cell or both
    SAVE_RESULTS = Read_Params('SAVE_RESULTS');
    FNames = Read_FNames;
    % CaData = [];
end


% Linux vs. Windows
if(ispc)
    slash = '\';
else
    slash = '/';
end



%% For each file in each folder, compute either whole image or individual
%% cell intensity depending on flag set by ANALYSIS_TYPE
Data = cell(size(FNames,1),1);
for i = 1:length(FNames) % folders
    Folder = FNames{i};
    if(iscell(Folder))
        Folder = Folder{1};
    end
    
    if(nargin<2)
    files = dir([SOURCE_PATH Folder slash '*.tif']);
    end
    % Remove files that are part of the exlusion list and remove
    % those whose name begins with Segmentation
    if(nargin==0)
        Exclude = Parse_FNames(FNames(i,:),'Exclude');
        ExcludeIdx = [];
        for j=1:numel(files)
            
            if(sum(strncmp(files(j).name,Exclude,3)))
                ExcludeIdx = [ExcludeIdx j];
            end
        end
        files(ExcludeIdx) = [];
    else
        exclude_idx = [];
        for j=1:numel(files)
            if(~isempty(strfind(files(j).name,'Segmentation')) || ...
                    ~isempty(strfind(files(j).name, 'MAX')))
                exclude_idx = [exclude_idx j];
            end
        end
        files(exclude_idx) = [];
    end
    
    % Make sure a segmentation file exists for each file. Don't need a
    % segmentation file for continuation images
    msg = '';
    disp_err=0;
    for j=1:numel(files)
        if(isempty(strfind(files(j).name,'-file')))
            if(~exist([SOURCE_PATH Folder slash 'Segmentation-' files(j).name(1:end-4) '.mat'],'file'))
                msg = sprintf('%s\n%s\n',msg,[SOURCE_PATH Folder slash 'Segmentation-' files(j).name]);
                disp_err=1;
            end
        end
    end
    if(disp_err)
        errordlg(msg,'Segmentation missing for these files');
        return
    end
    % Now we are ready to process these files
    % First make sure all files are readable.
    
    CaData = [];
    if(strcmp(ANALYSIS_TYPE,'BOTH'))
        L = cell(numel(files),1);
        info = cell(numel(files),1);
        for j=1:numel(files)
            
            filename = [SOURCE_PATH Folder slash files(j).name];
            try
                disp(['Reading file info for stack ' filename]);
                info{j} = imfinfo(filename);
            catch
                disp(['COULD NOT READ ' filename '. Cancel this function Ctrl+C and rerun']);
                pause
            end
            frames(j) = size(info{j},1);
            % Large files are split into filename-file001.tif,
            % -file002.tif etc
            idx = strfind(files(j).name,'file');
            if(idx)
                % This file is a continuation of a larger file. Use the
                % segmentation of the first file.
                stackname = files(j).name;
                stackname = [stackname(1:idx-2) '.tif'];
            else
                stackname = files(j).name;
            end
            % ** Segmentation file for a stack named "filename.tif" must be
            % ** named "Segmentation-filename.tif"
            segname = [SOURCE_PATH Folder slash 'Segmentation-' stackname(1:end-4) '.mat'];
            try
                disp(['Reading ' segname]);
                %                 RemoveNeurites(segname);
                data=load(segname);
                L{j} = data.L;
                if(max(max(L{j})) == 1)
                    L{j} = bwlabel(L{j});
                end
            catch
                disp(['COULD NOT READ ' segname '. Cancel this function Ctrl+C and rerun']);
                pause
            end
        end
    end
    
    % Create waitbar
    multiWaitbar('CloseAll','Name','');
    for j=1:numel(files)
        multiWaitbar(files(j).name,0,'Name','Gathering ROI-based fluorescence');
    end
    
    
    for j=1:numel(files) % filenames
        
        filename = [SOURCE_PATH Folder slash files(j).name];
        
        ImgInfo = info{j};
        switch ANALYSIS_TYPE
            case 'ExtractFirstFrame'
                disp(['Extracting first frame for ' filename]);
                I = imread(filename);
                imwrite(I,[strtok(filename,'.') '_first.tif'],'tif','Compression','none');
            case 'WHOLE'
                
                disp(['Computing full frame intensity for ' filename]);
                
                F = zeros(frames(j),1); % Store average image intensity
                
                for k=1:frames(j)
                    I = imread(filename,k);
                    F(k) = mean(mean(I));
                end
                
                CaData(j).filename = filename;
                CaData(j).F_whole = F;
              % Load fps from params.mat
              try
                  load('params.mat');
              catch
                  errordlg('params.mat does not exist. Go to Analyze -> Preferences -> Revert to default or set preferences','File does not exist','modal');
                  return;
              end
              CaData(j).fps = params.fps;
               
            case 'CELL'
                
                N = max(max(L{j})); % Number of cells
                F = zeros(N,frames(j));
                disp(['Computing cell based intensity for ' filename]);
                for k=1:frames(j)
                    I = imread(filename,k);
                    stats = regionprops(L{j},I,'MeanIntensity');
                    F(:,k) = [stats.MeanIntensity]';
                    
                end
                
                CaData(j).filename = filename;
                CaData(j).F_cell = F;
                CaData(j).L = L{j}
             % Load fps from params.mat
              try
                  load('params.mat');
              catch
                  errordlg('params.mat does not exist. Go to Analyze -> Preferences -> Revert to default or set preferences','File does not exist','modal');
                  return;
              end
              CaData(j).fps = params.fps;
            case 'BOTH'
                
                disp(['Processing ' filename ' for both avg image intensity and individual cell intensities']);
                
                N = max(max(L{j}));
                F_cell = zeros(N,frames(j));
                F_whole = zeros(frames(j),1);
                
                Icomposite = zeros(info{j}(1).Height,info{j}(1).Width);
                
                for k=1:frames(j)
                    if(~mod(k,100))
                        multiWaitbar(files(j).name,k/frames(j));
                    end
                    I = imread(filename,k); I = double(I);
                    F_whole(k) = mean(mean(I));
                    stats = regionprops(L{j},I,'MeanIntensity');
                    F_cell(:,k) = [stats.MeanIntensity]';
                    
                    Icomposite = Icomposite + I;
                end
                Icomposite = Icomposite/frames(j);
                image = zeros([size(Icomposite) 3],'uint16');
                image(:,:,2) = imadjust(uint16(Icomposite));
                
                CaData(j).image = image;
                CaData(j).filename = filename;
                CaData(j).F_cell = F_cell;
                CaData(j).F_whole = F_whole;
                CaData(j).L = L{j};
                multiWaitbar(files(j).name,1,'Color','g');
                % Load fps from params.mat
              try
                  load('params.mat');
              catch
                  errordlg('params.mat does not exist. Go to Analyze -> Preferences -> Revert to default or set preferences','File does not exist','modal');
                  return;
              end
              CaData(j).fps = params.fps;
                
        end
        analysis = CaData(j);
        [x,y,~] = fileparts(filename);
        save([x '/CaSignal-' y '.mat'],'analysis');
    end
    Data{i} = CaData;
    % Save CaData in same folder as the images
    if(SAVE_RESULTS)
        savefile = [TARGET_PATH Folder slash 'analysis.mat'];
        disp(['Saving results to ' savefile]);
        analysis = CaData(1:numel(files));
        save(savefile, 'analysis');
    end
end
