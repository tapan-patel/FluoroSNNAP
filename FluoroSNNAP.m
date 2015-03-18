function varargout = FluoroSNNAP(varargin)
% FluoroSNNAP MATLAB code for FluoroSNNAP.fig
%      FluoroSNNAP, by itself, creates a new FluoroSNNAP or raises the existing
%      singleton*.
%
%      H = FluoroSNNAP returns the handle to a new FluoroSNNAP or the handle to
%      the existing singleton*.
%
%      FluoroSNNAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FluoroSNNAP.M with the given input arguments.
%
%      FluoroSNNAP('Property','Value',...) creates a new FluoroSNNAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FluoroSNNAP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FluoroSNNAP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FluoroSNNAP

% Last Modified by GUIDE v2.5 15-Mar-2015 19:28:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FluoroSNNAP_OpeningFcn, ...
    'gui_OutputFcn',  @FluoroSNNAP_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FluoroSNNAP is made visible.
function FluoroSNNAP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FluoroSNNAP (see VARARGIN)

% Choose default command line output for FluoroSNNAP
handles.output = hObject;
handles.folders = [];
handles.curr_folder_idx = 0;
handles.curr_file_idx = 0;
handles.flims = []; % frame limits

handles.fps = [];
handles.time = [];
handles.imfinfo = [];
handles.Istack = [];
% Add paths
addpath([pwd '/FluoroSNNAP_code'],[pwd '/FluoroSNNAP_code/Cellsort/Cellsort1.3'],[pwd '/FluoroSNNAP_code/oopsi-master'],[pwd '/FluoroSNNAP_code/te_matlab_0'],[pwd '/FluoroSNNAP_code/mvgc_v1.0']);
% Try a TE job; if it fails, the user does not have a mex file
try
    asdf{1} = randsample(15,4);
    asdf{2} = randsample(15,4);
    asdf(end+1) = {10};
    asdf(end+1) = {[2,15]};
    [maxte, ci] = ASDFTE(asdf, 1:3);
catch
    hmsg=msgbox('***transent mex file does not exist. Attempting to compile now.***');
    try
        mex('FluoroSNNAP_code/te_matlab_0/transent.c');
    catch
        try
            delete(hmsg);
        end
        hmsg2=msgbox('Mex does not appear to be configured. Please run "mex -setup" from Matlab command window. Until this problem is fixed, you will not be able to use "transfer entropy" as a method for inferring functional connectivity.');
        uiwait(hmsg2);
    end
    try
        delete(hmsg);
    end
end

% Ask user if parallel computing toolbox is to be used
choice = questdlg('Do you want to use parallel computing? Requires proper installation of the parallel computing toolbox.','Enable parallel computing?','Enable parallel computing','Disable parallel computing','Enable parallel computing');
load('params.mat');
switch choice
    case 'Enable parallel computing'
        params.parallel = 1;
    case 'Disable parallel computing'
        params.parallel = 0;
end
save('FluoroSNNAP_code/params.mat','params');
if(params.parallel)
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        poolsize = 0;
    else
        poolsize = poolobj.NumWorkers;
    end
    if(poolsize==0)
        parpool;
    end
end
granger_set_paths;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FluoroSNNAP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FluoroSNNAP_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
idx = get(hObject,'Value');
if(isempty(idx))
    errordlg('Please add folders before making a selection','Bad input','modal');
    return;
end
handles.curr_folder_idx = idx;
handles.curr_folder = handles.folders{idx};
% Get tiff files from this folder and populate listbox1
fnames = dir([handles.curr_folder '/*.tif']);
handles.files = cell(length(fnames),1);
cntr = 1;
for i=1:length(fnames)
    if(isempty(strfind(fnames(i).name,'Segmentation-')))
        handles.files{cntr} = fnames(i).name;
        cntr = cntr+1;
    end
end
handles.files = handles.files(1:cntr-1);
set(handles.listbox2,'String',handles.files);
% If there are filenames with the word "-file" in it, warn the user that
% experiment tiff stacks are in separate files and it is recommended to
% merge them
if(~isempty(cell2mat(strfind(handles.files,'-file'))))
    msgbox(sprintf('Detected filenames with the pattern "-file".\nThis may represent continuation of an experiment\nUse File -> Merge before proceeding.'));
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
idx = get(hObject,'Value');
if(isempty(idx))
    errordlg('Bad selection. Please try again','Bad input','modal');
    return;
end
handles.curr_file_idx = idx;
handles.curr_file = handles.files{idx};

handles.fps = [];
handles.time = [];
handles.imfinfo = [];
handles.Istack = [];
handles.flims = [];
handles.Imean = [];
handles.Imax = [];
handles.data = [];

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1. Add folders
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt for folder
folder = uigetdir(pwd,'Select a folder containing tiff stacks');
handles.folders{end+1,1} = folder;
% Update listbox1 with the folder names
set(handles.listbox1,'String',handles.folders);
guidata(hObject,handles);
% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
% New - deletes current folders and files
function new_Callback(hObject, eventdata, handles)
% hObject    handle to new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folders = [];
handles.files = [];
handles.fps = [];
handles.time = [];
handles.imfinfo = [];
handles.Istack = [];
handles.flims = [];
set(handles.listbox1,'String',handles.folders);
set(handles.listbox2,'String',handles.files);
guidata(hObject,handles);

% --------------------------------------------------------------------
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton1_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function preprocessing_Callback(hObject, eventdata, handles)
% hObject    handle to preprocessing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function viewstack_Callback(hObject, eventdata, handles)
% hObject    handle to viewstack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    
    return;
end

% Ask for frame range
if(isempty(handles.imfinfo))
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    info = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    handles.imfinfo = info;
    multiWaitbar('Close All','Name','');
end
if(isempty(handles.flims))
    try
        answer = inputdlg({'First frame','Last frame'},'Select frame range',1,{'1',num2str(length(info))});
        
        handles.flims = [str2num(answer{1}),str2num(answer{2})];
    catch
        return;
    end
end
if(handles.flims(2)>length(handles.imfinfo))
    errordlg('Last frame cannot exceed the total number of frames in the image stack.','Bad Input','modal');
    
    return;
end

if(isempty(handles.Istack))
    multiWaitbar('CloseAll','Name','');
    multiWaitbar(['Loading ' handles.curr_file],'Busy','Name','Please wait');
    
    [Istack,time,fps] = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file),handles.flims);
    handles.Istack = Istack;
    handles.fps = fps;
    handles.time = time;
    try
        multiWaitbar('CloseAll','Name','');
    end
end
ViewImageStack(handles.Istack,handles.time,handles.fps);
guidata(hObject,handles);
% --------------------------------------------------------------------
function segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    
    return;
end
choice = questdlg('Do you want to use time-averaged image for segmentation (recommended) or use a single frame','Select a frame','Time-averaged image','select a frame','Time-averaged image');
switch choice
    case 'Time-averaged image'
        try
            if(isfield(handles,'Imean') && ~isempty(handles.Imean))
                ;
            else
                if(isempty(handles.imfinfo))
                    multiWaitbar('Close All','Name','');
                    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
                    info = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
                    handles.imfinfo = info;
                    multiWaitbar('Close All','Name','');
                end
                if(isempty(handles.flims))
                    try
                        answer = inputdlg({'First frame','Last frame'},'Select frame range',1,{'1',num2str(length(handles.imfinfo))});
                        
                        handles.flims = [str2num(answer{1}),str2num(answer{2})];
                    catch
                        return;
                    end
                end
                if(handles.flims(2)>length(handles.imfinfo))
                    errordlg('Last frame cannot exceed the total number of frames in the image stack.','Bad Input','modal');
                    
                    return;
                end
                
                multiWaitbar('Close All','Name','');
                multiWaitbar(['Loading ' handles.curr_file],0,'Name','Please wait');
                
                I = imread(fullfile(handles.curr_folder,handles.curr_file),...
                    'Info',handles.imfinfo,'Index',handles.flims(1));
                for i=handles.flims(1)+1:handles.flims(2)
                    if(~mod(i,100))
                        multiWaitbar(['Loading ' handles.curr_file],i/handles.flims(2),'Name','Please wait');
                    end
                    I = double(I) + double(imread(fullfile(handles.curr_folder,handles.curr_file),...
                        'Info',handles.imfinfo,'Index',i));
                end
                I = I./numel(handles.flims(1):handles.flims(2));;
                handles.Imean = I;
                multiWaitbar('Close All','Name','');
            end
        catch
            errordlg('Make a file selection first.','Bad input','modal');
            return
            
        end
    case 'select a frame'
        answer = inputdlg('Enter a frame number to use for segmentation (use ViewImageStack to find the frame #)','Select a frame',1,{'1'});
        idx = str2num(answer{1});
        try
            handles.Imean = imread(fullfile(handles.curr_folder,handles.curr_file),idx);
        catch
            errordlg('Make a file selection first.','Bad input','modal');
            return
        end
end

% % Ask for frame range
% if(isempty(handles.imfinfo))
%     info = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
%     handles.imfinfo = info;
% end
% if(isempty(handles.flims))
%     try
%         answer = inputdlg({'First frame','Last frame'},'Select frame range',1,{'1',num2str(length(info))});
%
%         handles.flims = [str2num(answer{1}),str2num(answer{2})];
%     catch
%         return;
%     end
% end
% if(handles.flims(2)>length(handles.imfinfo))
%     errordlg('Last frame cannot exceed the total number of frames in the image stack.','Bad Input','modal');
%
%     return;
% end
%
% if(isempty(handles.Istack))
%     multiWaitbar('CloseAll','Name','');
%     multiWaitbar(['Loading ' handles.curr_file],'Busy','Name','Please wait');
%
%     [Istack,time,fps] = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file),handles.flims);
%     handles.Istack = Istack;
%     handles.fps = fps;
%     handles.time = time;
%     try
%         multiWaitbar('CloseAll','Name','');
%     end
% end

SegmentationGUI(handles.Imean,fullfile(handles.curr_folder,handles.curr_file));
guidata(hObject,handles);
% --------------------------------------------------------------------
function meanimage_Callback(hObject, eventdata, handles)
% hObject    handle to meanimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    
    return;
end
if(isempty(handles.imfinfo))
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    info = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    handles.imfinfo = info;
    multiWaitbar('Close All','Name','');
    
end
if(isempty(handles.flims))
    try
        answer = inputdlg({'First frame','Last frame'},'Select frame range',1,{'1',num2str(length(info))});
        
        handles.flims = [str2num(answer{1}),str2num(answer{2})];
    catch
        return;
    end
end
if(handles.flims(2)>length(handles.imfinfo))
    errordlg('Last frame cannot exceed the total number of frames in the image stack.','Bad Input','modal');
    
    return;
end

multiWaitbar('Close All','Name','');
multiWaitbar(['Loading ' handles.curr_file],0,'Name','Please wait');

I = imread(fullfile(handles.curr_folder,handles.curr_file),...
    'Info',handles.imfinfo,'Index',handles.flims(1));
for i=handles.flims(1)+1:handles.flims(2)
    if(~mod(i,100))
        multiWaitbar(['Loading ' handles.curr_file],i/handles.flims(2),'Name','Please wait');
    end
    I = double(I) + double(imread(fullfile(handles.curr_folder,handles.curr_file),...
        'Info',handles.imfinfo,'Index',i));
end
I = I./numel(handles.flims(1):handles.flims(2));
handles.Imean = I;
multiWaitbar('Close All','Name','');

% if(isempty(handles.Istack))
%     multiWaitbar('CloseAll','Name','');
%     multiWaitbar(['Loading ' handles.curr_file],'Busy','Name','Please wait');
%
%     [Istack,time,fps] = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file),handles.flims);
%     handles.Istack = Istack;
%     handles.fps = fps;
%     handles.time = time;
%     try
%         multiWaitbar('CloseAll','Name','');
%     end
% end
% handles.Imean = mean(handles.Istack,3);
figure;imshow(imadjust(uint16(handles.Imean)),[]);
title(['Time-averaged image for ' fullfile(handles.curr_folder,handles.curr_file)]);
guidata(hObject,handles);
% --------------------------------------------------------------------
function maxprojection_Callback(hObject, eventdata, handles)
% hObject    handle to maxprojection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    
    return;
end

if(isempty(handles.imfinfo))
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    info = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    handles.imfinfo = info;
    multiWaitbar('Close All','Name','');
    
end
if(isempty(handles.flims))
    try
        answer = inputdlg({'First frame','Last frame'},'Select frame range',1,{'1',num2str(length(info))});
        
        handles.flims = [str2num(answer{1}),str2num(answer{2})];
    catch
        return;
    end
end
if(handles.flims(2)>length(handles.imfinfo))
    errordlg('Last frame cannot exceed the total number of frames in the image stack.','Bad Input','modal');
    
    return;
end
multiWaitbar('Close All','Name','');
multiWaitbar(['Loading ' handles.curr_file],0,'Name','Please wait');

I = imread(fullfile(handles.curr_folder,handles.curr_file),...
    'Info',handles.imfinfo,'Index',handles.flims(1));
for i=handles.flims(1)+1:handles.flims(2)
    if(~mod(i,100))
        multiWaitbar(['Loading ' handles.curr_file],i/handles.flims(2),'Name','Please wait');
    end
    Inew = imread(fullfile(handles.curr_folder,handles.curr_file),...
        'Info',handles.imfinfo,'Index',i);
    I(Inew>I) = Inew(Inew>I);
end

handles.Imax = I;
multiWaitbar('Close All','Name','');
% if(isempty(handles.Istack))
%     multiWaitbar('CloseAll','Name','');
%     multiWaitbar(['Loading ' handles.curr_file],'Busy','Name','Please wait');
%
%     [Istack,time,fps] = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file),handles.flims);
%     handles.Istack = Istack;
%     handles.fps = fps;
%     handles.time = time;
%     try
%         multiWaitbar('CloseAll','Name','');
%     end
% end
% handles.Imax = max(handles.Istack,[],3);
figure; imshow(imadjust(handles.Imax),[]);
title(['Maximum projection image for ' fullfile(handles.curr_folder,handles.curr_file)]);

guidata(hObject,handles);
% --------------------------------------------------------------------
function batchsegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to batchsegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Select a method for automated segmentation. ICA-based method is memory intensive and will take a long time (but probably more accurate). Active contour method is quick but will require manual touch-ups','Segmentation method','Active contour','ICA','Active contour');
switch choice
    case 'ICA'
        try
            prompt = {'Frames to process:','Number of principle components:', 'Number of PCs to use:','Weight of temporal information (0=pure spatial ICA, 1=pure temporal ICA)',...
                'Number of independent components','Standard deviation of Gaussian smoothing kernel (px)','Threshold for spatial filters (s.d.)'...
                'Min and max area of ROIs (px)','Downsample, scalar [temporal spatial]:','Visualize processing (slows computation)'};
            dlg_title = 'Input parameters for ICA-based segmentation';
            num_lines = 1;
            
            def = {'[1 1200]','150','1:100','1','100','4','2','[50 500]','[1 1]','0'};
            
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            
            flims = str2num(answer{1});
            nPCs = str2num(answer{2});
            PCuse = str2num(answer{3});
            mu = str2num(answer{4});
            nIC = str2num(answer{5});
            smwidth = str2num(answer{6});
            thresh = str2num(answer{7});
            arealim = str2num(answer{8});
            dsamp = str2num(answer{9});
            plotting = str2num(answer{10});
            
            badframes = [];
            outputdir = [];
        catch
            return;
        end
        dsamp = [];
        badframes = [];
        
        for i=1:length(handles.folders)
            fnames = dir([handles.folders{i} '/*.tif']);
            handles.files = [];
            cntr = 1;
            for j=1:length(fnames)
                if(isempty(strfind(fnames(j).name,'Segmentation-')))
                    handles.files{cntr} = fnames(j).name;
                    cntr = cntr+1;
                end
            end
            folder = handles.folders{i};
            files = handles.files;
            for j=1:length(handles.files)
                % Batch each of these files for ICA-based segmentation
                fn = fullfile(folder,files{j});
                disp(['Processing ' fn]);
                
                info = imfinfo(fn);
                
                ICASegmentation(fn,[1 length(info)],nPCs,PCuse,mu,nIC,smwidth,thresh,arealim,plotting,dsamp,badframes,handles.folders{i});
            end
        end
        multiWaitbar('CloseAll','Name','');
    case 'Active contour'
        for i=1:length(handles.folders)
            fnames = dir([handles.folders{i} '/*.tif']);
            handles.files = [];
            cntr = 1;
            for j=1:length(fnames)
                if(isempty(strfind(fnames(j).name,'Segmentation-')))
                    handles.files{cntr} = fnames(j).name;
                    cntr = cntr+1;
                end
            end
            folder = handles.folders{i};
            files = handles.files;
            multiWaitbar('Close all','Name','');
            for j=1:length(handles.files)
                multiWaitbar(handles.files{j},0,'Name','Active contour segmentation');
            end
            for j=1:length(handles.files)
                
                % Batch each of these files for active contour segmentation
                fn = fullfile(folder,files{j});
                info = imfinfo(fn);
                frames = numel(info);
                Imean = double(imread(fn,1));
                for l=2:frames
                    if(~mod(l,100))
                        multiWaitbar(handles.files{j},l/frames*.9);
                    end
                    I = imread(fn,l);
                    Imean = Imean + double(I);
                end
                Imean = Imean./frames;
                BW = Imean>prctile(Imean(:),95);
                BW = activecontour(Imean,BW);
                L = bwlabel(BW);
                C = regionprops(L,'Area');
                A = [C.Area];
                for l=1:numel(A)
                    if(A(l)<10)
                        L(L==l)=0;
                    end
                end
                L = bwlabel(L);
                ica=0;
                
                savename = [folder '/Segmentation-' files{j}(1:end-4) '.mat'];
                save(savename,'L','ica');
                multiWaitbar(handles.files{j},1,'Color','g');
            end
            multiWaitbar('Close all','Name','');
        end
end
% --------------------------------------------------------------------
function viewsegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to viewsegmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    %      uicontrol(handles.figure1)
    return;
end
try
    [~,prefix] = fileparts(handles.curr_file);
    fname = [handles.curr_folder '/Segmentation-' prefix '.mat'];
    data = load(fname);
    figure;
    imshow(label2rgb(data.L));
    
    axis image;
    title([num2str(max(data.L(:))) ' total ROIs']);
    hold on
    C = regionprops(data.L,'Centroid');
    for i=1:length(C)
        text(C(i).Centroid(1),C(i).Centroid(2),num2str(i),'Color','k');
    end
    figure;
    imagesc(data.L); title('Labeled segmentation'); colorbar; axis image;
    colormap('jet');
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
catch
    errordlg(['Could not load ' fname]);
    return
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
batchsegmentation_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function analysis_Callback(hObject, eventdata, handles)
% hObject    handle to analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function batchanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to batchanalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.folders = get(handles.listbox1,'String');

AnalyzeData(handles.folders);
PostProcess(handles.folders);



% --------------------------------------------------------------------
function single_file_Callback(hObject, eventdata, handles)
% hObject    handle to single_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

AnalyzeData({handles.curr_folder},handles.curr_file);
PostProcess({handles.curr_folder},handles.curr_file);



% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function caltrace_Callback(hObject, eventdata, handles)
% hObject    handle to caltrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Plot the fluorescence trace for all neurons on a single plot
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
figure;
t = 0:1/data.fps:size(data.F_cell,2)/data.fps-1/data.fps;
sig = data.F_cell;

% calculate shift
mi = min(sig,[],2);
ma = max(sig,[],2);
shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,1,size(sig,2));

%plot data
plot(t,sig+shift)

% edit axes
set(gca,'ytick',mean(sig+shift,2),'yticklabel',1:data.N)
grid on
ylim([mi(1) max(max(shift+sig))]);
xlabel('Time (s)'); ylabel('Neuron ID');
title(data.filename);
% --------------------------------------------------------------------
function raster_Callback(hObject, eventdata, handles)
% hObject    handle to raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end

figure;
hold on
for i=1:data.N
    spks = data.Spikes_cell{i};
    if(~isempty(spks))
        plot(spks/data.fps,i,'b.');
    end
end
xlabel('Time (s)'); ylabel('Neuron ID');
title(['Calcium activity for ' data.filename]);
hold off
% --------------------------------------------------------------------
function spikeinspect_Callback(hObject, eventdata, handles)
% hObject    handle to spikeinspect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');

EventInspector(data);

% --------------------------------------------------------------------
function FC_Callback(hObject, eventdata, handles)
% hObject    handle to FC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end


choice = questdlg('Do you want to use pre-processed data or recompute?','Choice','Pre-processed','Recompute','Pre-processed');
switch choice
    case 'Recompute'
        y=FC_params;
        hmsg = msgbox('Please select a method for inferring functional connectivity and update parameters as needed');
        uiwait(hmsg);
        uiwait(y);
        recompute_Callback(hObject,[],handles);
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
s = handles.data;

params = s.params;
try
    L = s.L;
    N = s.N;
    C = s.SynchroCluster;
    if(params.FC.method_idx==1)
        A = s.FC.CC.A;
        Ci = s.FC.CC.modularity_Ci;
        Q = s.FC.CC.modularity_Q;
        clu = s.FC.CC.clustering_coef;
    elseif(params.FC.method_idx==2)
        A = s.FC.PC.A;
        Ci = s.FC.PC.modularity_Ci;
        Q = s.FC.PC.modularity_Q;
        clu = s.FC.PC.clustering_coef;
    elseif(params.FC.method_idx==3)
        A = s.FC.phase.A;
        Ci = s.FC.phase.modularity_Ci;
        Q = s.FC.phase.modularity_Q;
        clu = s.FC.phase.clustering_coef;
    elseif(params.FC.method_idx==4)
        A = s.FC.GC.A;
        Ci = s.FC.GC.modularity_Ci;
        Q = s.FC.GC.modularity_Q;
        clu = s.FC.GC.clustering_coef;
    elseif(params.FC.method_idx==5)
        A = s.FC.TE.peakTE;
        Ci = s.FC.TE.modularity_Ci;
        Q = s.FC.TE.modularity_Q;
        clu = s.FC.TE.clustering_coef;
    else
        A = s.FC.phase.A;
        Ci = s.FC.phase.modularity_Ci;
        Q = s.FC.phase.modularity_Q;
        clu = s.FC.phase.clustering_coef;
    end
    
    c = regionprops(L,'Centroid');
    c = reshape([c.Centroid],2,N);
    x = c(1,:);
    y = c(2,:);
    DF = s.DF; rise_time = s.rise_time; fall_time = s.fall_time; CV = s.CV;
    
    SI = s.SI_m;
    Colors = varycolor(max(Ci));
catch
    errordlg('Could not use pre-processed information. Missing entries. Select "Reprocess" option');
    return;
end


netfig = figure;

hdt = datacursormode(netfig);
set(hdt,'DisplayStyle','window');

set (netfig, 'Units', 'normalized', 'Position', [0,0,1,1]);
h1 = subplot(2,2,1);
set(hdt,'UpdateFcn',{@labeltip,C,A,c,s.Spikes_cell,DF,...
    rise_time,fall_time,clu,h1});
if(size(s.image,3)~=3)
    image = zeros([size(s.image) 3]);
    image(:,:,2) = imadjust(s.image);
else
    image = s.image;
end
imagesc(image,'Parent',h1);
hold on

gplot(A,[x' y']);

set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);

hold on
plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
set(gca,'YDir','reverse');
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
title({['N = ' num2str(N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});


subplot(2,2,2);
[~,pos] = max(C.PI,[],2);
[~,IDX] = sort(pos);
imagesc(C.C(IDX,IDX),[0 1]); axis square; colorbar;
colormap('jet');
xlabel('Neuron ID (rearranged)'); ylabel('Neuron ID (rearranged)');
msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
title(msg);

subplot(2,2,3);

set(gca, 'ColorOrder', Colors);
hold all;
for i=1:max(Ci)
    scatter(x(Ci==i),y(Ci==i),'filled','SizeData',10^2);
end
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca,'YDir','reverse');

%                 title({sprintf('N = %d. Modules = %d, Modularity = %f',length(x),max(Ci),Q)});

title_string = sprintf('N = %d. Modules = %d. Modularity = %f\n',length(x),max(Ci),Q);
title_string = sprintf('%sSynchronization by modules: ',title_string);
for k=1:length(SI)
    title_string = sprintf('%s %f,',title_string,SI(k));
end
title(title_string);
subplot(2,2,4);
t = 0:1/s.fps:size(s.F_cell,2)/s.fps-1/s.fps;
hold all;
cntr = 1;
for i=1:max(Ci)
    modules = find(Ci==i);
    for j=1:length(modules)
        try
            spks = s.Spikes_cell{modules(j)};
            plot(spks/s.fps,cntr,'.','Color',Colors(i,:));
            cntr = cntr+1;
        end
    end
end
title('Activity rearranged by modules');


function out_text = labeltip(h,e,C,A,c,S,DF,rise_time,fall_time,clu,h1)
N = size(A,1);

pos = get(e,'Position');
neuron = find(c(1,:)==pos(1) & c(2,:)==pos(2));
if(isempty(neuron))
    out_text = '';
    plot(c(1,:),c(2,:),'ro','MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
    plot(c(1,~logical(sum(A)) & ~logical(sum(A'))),c(2,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
else
    
    % Reset the size and color of all neurons to the way it was
    % originally
    
    %                     cla(h1);
    %                    imagesc(image,'Parent',h1);
    %                 hold on
    %
    %                 gplot(A,[x' y']);
    %
    %                 set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
    %                 set(gca,'XTickLabel',[]);
    %                 set(gca,'YTickLabel',[]);
    %                 hold on
    %                 plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
    %                 set(gca,'YDir','reverse');
    %                 set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    %                 title({['N = ' num2str(N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});
    
    
    out_text = ['Neuron #: ' num2str(neuron)];
    out_text = sprintf('%s. CI: %.02f',out_text,sum(A(neuron,:))/N);
    out_text = sprintf('%s,   PI: %.02f,  CC: %.02f',out_text,max(C.PI(neuron,:)),clu(neuron));
    
    out_text = sprintf('%s. Events: %d, Amplitude: %.02f, <t_rise> = %.02f (s), <t_fall> = %.02f (s)',...
        out_text,length(S{neuron}),DF(neuron),mean(rise_time{neuron}),...
        mean(fall_time{neuron}));
    
    connecting = find(A(neuron,:));
    hold on
    plot(c(1,:),c(2,:),'ro','MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);
    plot(c(1,~logical(sum(A)) & ~logical(sum(A'))),c(2,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
    plot(c(1,connecting),c(2,connecting),'bo','MarkerFaceColor','b','MarkerSize',6)
    
end

% --------------------------------------------------------------------
function synchronization_Callback(hObject, eventdata, handles)
% hObject    handle to synchronization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

msg1 = sprintf('Synchronization cluster analysis method\n');
msg1 = sprintf('%s1 = phase, 2 = equal-time correlation, 3 = entropy',msg1);
msg2 = sprintf('Number of times to perform surrogate resampling:');
msg3 = sprintf('Minimum size of synchronization cluster (# of neurons): ');
try
    load params.mat
end
try
    answer = inputdlg({msg1,msg2,msg3},'Synchronization cluster analysis',1,{num2str(params.sca_type),num2str(params.sca_N),num2str(params.sca_size)});
catch
    answer = inputdlg({msg1,msg2,msg3},'Synchronization cluster analysis',1,{'1','20','3'});
end
params.sca_type = str2num(answer{1});
params.sca_N = str2num(answer{2});
params.sca_size = str2num(answer{3});
save('FluoroSNNAP_code/params.mat','params');
% --------------------------------------------------------------------
function grouped_raster_Callback(hObject, eventdata, handles)
% hObject    handle to grouped_raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
if(data.params.FC)
figure
hold all;
cntr = 1;
Ci = data.modules;
Colors = varycolor(max(Ci));
for q=1:max(Ci)
    modules = find(Ci==q);
    for m=1:length(modules)
        try
            evnts = data.Spikes_cell{modules(m)};
            plot(evnts/data.fps,cntr,'.','Color',Colors(q,:));
            cntr = cntr+1;
        end
    end
end
title('Activity rearranged by modules');
ylim([0 data.N]); xlim([0 size(data.F_cell,2)/data.fps]);
xlabel('Time (s)'); ylabel('Neuron ID (rearranged)');
else
    errordlg('You chose not to infer functional connectivity. Therefore, this option is not available to you. Please go to Analysis -> Preferences -> Choose analysis modules','Bad input','modal');
    return;
end


% --------------------------------------------------------------------
function coactive_Callback(hObject, eventdata, handles)
% hObject    handle to coactive (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
spk = data.foopsi;

% Binarize spks. Prompt user for number of s.d above 0
try
    answer = inputdlg({'Number of standard deviations above 0','% of co-active cells that make an ensemble'},'spike probability',1,{'3','20'});
    
    thresh = str2num(answer{1});
    frac_coactive = str2num(answer{2})/100;
catch
    return;
end
APs = zeros(size(spk));
for i=1:size(spk,1)
    sd = std(spk(i,:));
    APs(i,:) = spk(i,:)>sd*thresh;
end
simult = sum(APs,1);
t = 0:1/data.fps:size(data.F_cell,2)/data.fps-1/data.fps;
figure
plot(t,simult./data.N*100,'k','LineWidth',2); xlabel('Time (s)'); ylabel('% of cells co-active');
figure
hold on
for i=1:size(APs,2)
    idx = find(APs(:,i));
    if(~isempty(idx))
        if(length(idx)/data.N>frac_coactive)
            plot(t(i),idx,'r.');
        else
            plot(t(i),idx,'b.');
        end
    end
end
xlabel('Time (s)'); ylabel('Neuron ID'); title('Inferred activity. Co-active cells labeled in red');

% --------------------------------------------------------------------
function viewspikes_Callback(hObject, eventdata, handles)
% hObject    handle to viewspikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
if(data.params.analyze.spike_probability)
spk = data.foopsi;


figure;
t = 0:1/data.fps:size(data.F_cell,2)/data.fps-1/data.fps;
sig = spk;

% calculate shift
mi = min(sig,[],2);
ma = max(sig,[],2);
shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,1,size(sig,2));

%plot data
plot(t,sig+shift)

% edit axes
set(gca,'ytick',mean(sig+shift,2),'yticklabel',1:data.N)
grid on
ylim([mi(1) max(max(shift+sig))]);
xlabel('Time (s)'); ylabel('Neuron ID'); title('Inferred spikes');
else
    errordlg('You chose not to infer functional connectivity. Therefore, this option is not available to you. Please go to Analysis -> Preferences -> Choose analysis modules','Bad input','modal');
    return;
end

guidata(hObject,handles);
% --------------------------------------------------------------------
function inferred_raster_Callback(hObject, eventdata, handles)
% hObject    handle to inferred_raster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
if(data.params.analyze.spike_probability)
spk = data.foopsi;

% Binarize spks. Prompt user for number of s.d above 0
try
    answer = inputdlg('Number of standard deviations above 0','Threshold spike probability',1,{'3'});
    
    thresh = str2num(answer{1});
catch
    return;
end
APs = zeros(size(spk));
for i=1:size(spk,1)
    sd = std(spk(i,:));
    APs(i,:) = spk(i,:)>(sd*thresh);
end

% Plot binarized activity map
figure;
hold on
for i=1:data.N
    spks = find(APs(i,:));
    if(~isempty(spks))
        plot(spks/data.fps,i,'b.');
    end
end
xlabel('Time (s)'); ylabel('Neuron ID');
title(['Inferred spike activity for ' data.filename]);
hold off
else
    errordlg('You chose not to infer functional connectivity. Therefore, this option is not available to you. Please go to Analysis -> Preferences -> Choose analysis modules','Bad input','modal');
    return;
end


% --------------------------------------------------------------------
function foopsi_Callback(hObject, eventdata, handles)
% hObject    handle to foopsi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function preferences_analysis_Callback(hObject, eventdata, handles)
% hObject    handle to preferences_analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function cell_summary_Callback(hObject, eventdata, handles)
% hObject    handle to cell_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
SingleCellReport(data);


% --------------------------------------------------------------------
function template_library_Callback(hObject, eventdata, handles)
% hObject    handle to template_library (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function h=view_waveforms_Callback(hObject, eventdata, handles)
% hObject    handle to view_waveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    load('spikes.mat');
catch
    msg = sprintf('Could not load spikes.mat. Please make sure file exists\n');
    errordlg(msg,'Cannot locate file','modal');
    return;
end

n = length(spikes);
r = ceil(sqrt(n));
c = ceil(n/r);
h=figure;
cntr = 1;
for i=1:n
    subplot(r,c,cntr); plot(spikes{i},'k'); title(['Waveform #' num2str(i)]);
    cntr = cntr+1;
end
handles.waveforms = spikes;
guidata(hObject,handles);
% --------------------------------------------------------------------
function add_waveform_Callback(hObject, eventdata, handles)
% hObject    handle to add_waveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end
try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    msg = sprintf('Could not load %s. Please make sure file exists or re-run analysis\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end

% Ask for ROI to load
try
    answer = inputdlg('Enter ROI number to view its fluorescence trace','ROI number',1,{'1'});
    
    ROI = str2num(answer{1});
catch
    return;
end
h=figure;
plot(data.F_cell(ROI,:),'k');
xlabel('Frame #'); ylabel('Fluorescence (a.u.)');
title(sprintf('Select two points that mark the start and end of a waveform to add to the library\nPress <Enter> after selecting the second point'));
[x,~] = getpts(h);
load('spikes.mat');
spikes{end+1} = data.F_cell(ROI,floor(x(1)):floor(x(2)));
close(h);
figure; plot(spikes{end},'k');
msgbox('Waveform added to library');
save('spikes.mat','spikes');
guidata(hObject,handles);
% --------------------------------------------------------------------
function delete_waveform_Callback(hObject, eventdata, handles)
% hObject    handle to delete_waveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Show the waveforms
h=view_waveforms_Callback(hObject, eventdata, handles);
% Ask user for which waveforms they wish to remove
try
    answer = inputdlg('Enter waveforms you would like to delete','Delete waveforms',1,{'1,2'});
    
    idx = str2num(answer{1});
catch
    return;
end
handles.waveforms(idx) = [];
spikes = handles.waveforms;
save('spikes.mat','spikes');
% Show new waveforms
delete(h);
view_waveforms_Callback(hObject, eventdata, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Prompt for folder
folder = uigetdir(pwd,'Select a folder containing .csv files');
files = dir([folder '/*.csv']);
for i=1:numel(files)
    CSVtoTiff(fullfile(folder,files(i).name));
end
handles.folders{end+1,1} = folder;
% Update listbox1 with the folder names
set(handles.listbox1,'String',handles.folders);
guidata(hObject,handles);

% --------------------------------------------------------------------
function open_csv_Callback(hObject, eventdata, handles)
% hObject    handle to open_csv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pushbutton3_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function baseline_fluorescence_Callback(hObject, eventdata, handles)
% hObject    handle to baseline_fluorescence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg1 = sprintf('F0 = mean of the lower **%% of previous **-s values');
msg1 = sprintf('%s\nEnter the time window (seconds): ',msg1);
msg2 = sprintf('Enter percentile for baseline fluorescence: ');
try
    load params.mat
end
try
    answer = inputdlg({msg1,msg2},'F0',1,{num2str(params.F0_time),num2str(params.F0_pctl)});
catch
    answer = inputdlg({msg1,msg2},'F0',1,{'10','50'});
end
params.F0_time = str2num(answer{1});
params.F0_pctl = str2num(answer{2});
save('FluoroSNNAP_code/params.mat','params');

% --------------------------------------------------------------------
function event_detection_Callback(hObject, eventdata, handles)
% hObject    handle to event_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg1 = sprintf('Template-based = 1\nDeconvolution = 2');
msg2 = sprintf('Threshold: [0,1] for template-based; scalar >=1 for deconvolution: ');
msg3 = sprintf('Minimum deltaF/F amplitude: ');
try
    load params.mat
end
try
    answer = inputdlg({msg1,msg2,msg3},'Event detection',1,{num2str(params.event_type),num2str(params.event_thresh),num2str(params.event_amplitude)});
catch
    answer = inputdlg({msg1,msg2,msg3},'Event detection',1,{'1','0.85','0.01'});
end
params.event_type = str2num(answer{1});
params.event_thresh = str2num(answer{2});
params.event_amplitude = str2num(answer{3});
if(params.event_type==2 && params.event_thresh<1)
    errordlg('Threshold must be a scalar >=1 if deconvolution based event detection is selected','Bad input','modal');
    return
end
save('FluoroSNNAP_code/params.mat','params');


% --------------------------------------------------------------------
function functional_connectivity_Callback(hObject, eventdata, handles)
% hObject    handle to functional_connectivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

FC_params;
% --------------------------------------------------------------------
function defaults_Callback(hObject, eventdata, handles)
% hObject    handle to defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Revert to default parameters
params.parallel = 1;
params.fps = 10;
params.F0_time = 10; % 10 seconds
params.F0_pctl = 50;
params.event_type =1;
params.event_thresh = 0.85;
params.event_amplitude = 0.01;
params.alpha_level = 0.001;
params.sca_type = 1;
params.sca_N = 20;
params.sca_size = 3;
params.spike_inference_iter = 500;
params.spike_inference_filter = 1;
params.spike_inference_A = 50;
params.spike_inference_n = 1;
params.spike_inference_kd = 200;
params.FC.CC.Nsur = 100;
params.FC.CC.maxlag = .5;
params.FC.PC.alpha = 0.001;
params.FC.phase.Nsur = 100;
params.FC.phase.alpha = 0.001;
params.FC.GC.morder=20;
params.FC.GC.alpha=0.05;
params.FC.GC.iter=100;
params.FC.TE.lags=1:10;
params.FC.TE.Nsur=100;
params.FC.method_idx = 6;
params.network_ensemble_sd = 3;
params.network_ensemble_Nsur = 1000;
params.analyze.deltaF=1;
params.analyze.detect_events=1;
params.analyze.sca=1;
params.analyze.FC=1;
params.analyze.controllability=1;
params.analyze.kinetics=1;
params.analyze.spike_probability=1;
params.analyze.ensembles=1;
params.analyze.figure=1;

save('FluoroSNNAP_code/params.mat','params');

% --------------------------------------------------------------------
function acquisition_fps_Callback(hObject, eventdata, handles)
% hObject    handle to acquisition_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg = sprintf('Enter acquisition frame rate (Hz): ');
try
    load params.mat
end
try
    answer = inputdlg(msg,'FPS',1,{num2str(params.fps)});
catch
    answer = inputdlg(msg,'FPS',1,{'10'});
end
params.fps = str2num(answer{1});
save('FluoroSNNAP_code/params.mat','params');


% --------------------------------------------------------------------
function spike_inference_Callback(hObject, eventdata, handles)
% hObject    handle to spike_inference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg1 = sprintf('Max number of iterations');
msg2 = sprintf('Pre-process with high-pass filter (0=no, 1=yes)');
msg3 = sprintf('Jump in [Ca2+] transient after an action potential');
msg4 = sprintf('Hill coefficient:');
msg5 = sprintf('Kd (dissociation constant:');

try
    load params.mat
end
try
    try
        answer = inputdlg({msg1,msg2,msg3,msg4,msg5},'Spike inference',1,...
            {num2str(params.spike_inference_iter),num2str(params.spike_inference_filter),num2str(params.spike_inference_A),num2str(params.spike_inference_n),num2str(params.spike_inference_kd)});
    catch
        answer = inputdlg({msg1,msg2,msg3,msg4,msg5},'Spike inference',1,{'500','1','50','1','200'});
    end
    params.spike_inference_iter = str2num(answer{1});
    params.spike_inference_filter = str2num(answer{2});
    params.spike_inference_A = str2num(answer{3});
    params.spike_inference_n = str2num(answer{4});
    params.spike_inference_kd = str2num(answer{5});
    save('FluoroSNNAP_code/params.mat','params');
end
% --------------------------------------------------------------------
function merge_files_Callback(hObject, eventdata, handles)
% hObject    handle to merge_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt the user for indices, in order to merge
try
    handles.curr_folder_idx;
    handles.curr_folder;
catch
    errordlg('You must first select a folder containing tiff files','Bad input','modal');
    return;
end
msg = sprintf('Look at the listbox and enter the indices (in order) of the files you would like to merge.\nEnter the indices as [index1, index2, index3, etc]');
try
    answer = inputdlg(msg,'Enter indices',1,{'[2,1]'});
    
    indices = str2num(answer{1});
catch
    return;
end

for i=2:length(indices)
    hmsg = waitbar(0,sprintf('Please wait. Appending %s to %s\n.',handles.files{indices(i)},handles.files{indices(1)}));
    info = imfinf(fullfile(handles.curr_folder,handles.files{indices(i)}));
    
    frames = numel(info);
    for k=1:frames
        if(~mod(k,100))
            waitbar(k/frames,hmsg,sprintf('Please wait. Appending %s to %s\n.',handles.files{indices(i)},handles.files{indices(1)}));
        end
        I = imread(fullfile(handles.curr_folder,handles.files{indices(i)}),'Info',info,'Index',k);
        imwrite(I,fullfile(handles.curr_folder,handles.files{indices(1)}),'WriteMode','append');
    end
    try
        delete(hmsg);
    end
    % Recycle files that were appended
    recycle('on');
    delete(fullfile(handles.curr_folder,handles.files{indices(i)}));
    
end
% Reload the listbox
handles.files(indices(2:end))=[];
set(handles.listbox2,'String',handles.files);
guidata(hObject,handles);


% --------------------------------------------------------------------
function remove_frames_Callback(hObject, eventdata, handles)
% hObject    handle to remove_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
% Prompt user for frame #s to remove from image stack
try
    answer = inputdlg('Enter frame numbers that you would like to delete from the image stack','Frames to remove',1,{'200:300'});
    
    indices = str2num(answer{1});
catch
    return;
end
fname = fullfile(handles.curr_folder,handles.curr_file);
hmsg = waitbar(0,sprintf('Please wait while bad frames are removed and the image stack re-written\n%s',fname));
multiWaitbar('Close All','Name','');
multiWaitbar('Reading image stack info','Busy','Name','Please wait');
info = imfinfo(fname);
multiWaitbar('Close All','Name','');

frames = numel(info);
waitbar(.5,hmsg,sprintf('Please wait while bad frames are removed and the image stack re-written\n%s',fname));
recycle('on');
delete(fname);
for i=1:frames
    I = imread(fname,'Info',info,'Index',i);
    if(~nnz(i==indices))
        if(~mod(i,100))
            waitbar(.5+.5*i/size(Istack,3),hmsg,sprintf('Please wait while bad frames are removed and the image stack re-written\n%s',fname));
        end
        imwrite(I,fname,'WriteMode','append');
    end
end
try
    delete(hmsg);
end



% --------------------------------------------------------------------
function Untitled_4_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function check_artifact_Callback(hObject, eventdata, handles)
% hObject    handle to check_artifact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end

% Read the first and last frames and compute a 2D correlation coefficient.
% If r>0.7, it is likely the frames are aligned. Then pick 100 frames at
% random and compute r. If all r for all 100 is >0.7, fairly certain the
% frames are aligned
if(isfield(handles,'imfinfo') && ~isempty(handles.imfinfo))
    N = numel(handles.imfinfo);
else
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    handles.imfinfo = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('Close All','Name','');
    
    N = numel(handles.imfinfo);
end
I1 = imread(fullfile(handles.curr_folder,handles.curr_file),1);
I2 = imread(fullfile(handles.curr_folder,handles.curr_file),N);
r = corr2(I1,I2);
if(r<0.7)
    msg = sprintf('Suspect that %s contains motion artifact.\nRecommend scrolling through the stack using Pre-processing -> View Image Stack.',fullfile(handles.curr_folder,handles.curr_file));
    msgbox(msg,'Motion artifact suspected','warn');
    return
end

for i=1:100
    idx=randsample(N,2);
    I1 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(1));
    I2 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(2));
    r = corr2(I1,I2);
    if(r<0.7)
        msg = sprintf('Suspect that %s contains motion artifact.\nRecommend scrolling through the stack using Pre-processing -> View Image Stack.',fullfile(handles.curr_folder,handles.curr_file));
        msgbox(msg,'Motion artifact suspected','warn');
        return;
    end
end
msgbox('Image stack does not appear to have motion artifact. Run correlation matrix for a more detailed look');

guidata(hObject,handles);
% --------------------------------------------------------------------
function view_frames_Callback(hObject, eventdata, handles)
% hObject    handle to view_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
if(isfield(handles,'imfinfo') && ~isempty(handles.imfinfo))
    N = numel(handles.imfinfo);
else
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    handles.imfinfo = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('Close All','Name','');
    
    N = numel(handles.imfinfo);
end
try
    answer = inputdlg(sprintf('Enter two frames between 1 and %d to view overlay.',numel(handles.imfinfo)),'Enter two frames',1,{sprintf('[1 %d]',numel(handles.imfinfo))});
    idx = str2num(answer{1});
catch
    errordlg('Bad input format. Please try again.');
    return
end
I1 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(1));
I2 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(2));
r = corr2(I1,I2);
I1 = imadjust(I1);
I2 = imadjust(I2);
figure; imshowpair(I1,I2,'Scaling','joint'); title(sprintf('2-D correlation between frames %d and %d = %f',idx(1),idx(2),r));
guidata(hObject,handles);
% --------------------------------------------------------------------
function align_pair_frames_Callback(hObject, eventdata, handles)
% hObject    handle to align_pair_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
if(isfield(handles,'imfinfo') && ~isempty(handles.imfinfo))
    N = numel(handles.imfinfo);
else
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    handles.imfinfo = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('Close All','Name','');
    
    N = numel(handles.imfinfo);
end
try
    answer = inputdlg(sprintf('Enter two frames between 1 and %d to align.\nBy default, the first frame is the template',numel(handles.imfinfo)),'Enter two frames',1,{sprintf('[1 %d]',numel(handles.imfinfo))});
    idx = str2num(answer{1});
catch
    errordlg('Bad input format. Please try again.');
    return
end
[optimizer,metric] = imregconfig('monomodal');
registration.metric.MattesMutualInformation;
if(isfield(handles,'transformiter') && ~isempty(handles.transformiter))
    optimizer.MaximumIterations=handles.transformiter;
end
transform = 'translation';
if(isfield(handles,'transformtype') && ~isempty(handles.transformtype))
    switch transformtype
        case 1
            transform = 'traslation';
        case 2
            transform = 'rigid';
        case 3
            transform = 'similarity';
        case 4
            transform = 'affine';
    end
end
I1 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(1));
I2 = imread(fullfile(handles.curr_folder,handles.curr_file),idx(2));
I2_warped = imregister(I2,I1,transform,optimizer,metric);

r = corr2(I1,I2);
r_warped = corr2(I1,I2_warped);
% Show the orignal and warped image overlays
I1 = imadjust(I1);
I2 = imadjust(I2);
I2_warped = imadjust(I2_warped);
figure; subplot(1,2,1); imshowpair(I1,I2,'Scaling','joint'); axis image;
title(sprintf('ORIGINAL overlay between frames %d and %d.\n2-D correlation = %f',idx(1),idx(2),r));
subplot(1,2,2); imshowpair(I1,I2_warped,'Scaling','joint'); axis image;
title(sprintf('REGISTERED overlay between frames %d and %d.\n2-D correlation = %f',idx(1),idx(2),r_warped));
guidata(hObject,handles);
% --------------------------------------------------------------------
function align_to_template_Callback(hObject, eventdata, handles)
% hObject    handle to align_to_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
if(isfield(handles,'imfinfo') && ~isempty(handles.imfinfo))
    N = numel(handles.imfinfo);
else
    multiWaitbar('Close All','Name','');
    multiWaitbar('Reading image stack info','Busy','Name','Please wait');
    handles.imfinfo = imfinfo(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('Close All','Name','');
    
    N = numel(handles.imfinfo);
end
try
    answer = inputdlg({'Enter the template frame number',sprintf('Enter a range of frames between 1 and %d that will be aligned to the template frame',numel(handles.imfinfo))},'Enter frame range for registration',1,{'1','200:300'});
    template_idx = str2num(answer{1});
    idx = str2num(answer{2});
catch
    errordlg('Bad input format. Please try again.');
    return
end

hwait = waitbar(0,'Please wait while aligning image frames to template',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
% Read the image stack
if(isfield(handles,'Istack') && ~isempty(handles.Istack))
    ;
else
    multiWaitbar('CloseAll','Name','');
    multiWaitbar('Reading image stack.','Busy','Name','Please wait');
    handles.Istack = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('CloseAll','Name','');
end

Istack_warped = handles.Istack;
[optimizer,metric] = imregconfig('monomodal');
registration.metric.MattesMutualInformation;
if(isfield(handles,'transformiter') && ~isempty(handles.transformiter))
    optimizer.MaximumIterations=handles.transformiter;
end
transform = 'translation';
if(isfield(handles,'transformtype') && ~isempty(handles.transformtype))
    switch transformtype
        case 1
            transform = 'traslation';
        case 2
            transform = 'rigid';
        case 3
            transform = 'similarity';
        case 4
            transform = 'affine';
    end
end
for i=1:length(idx)
    if getappdata(hwait,'canceling')
        break
    end
    Istack_warped(:,:,idx(i)) = imregister(handles.Istack(:,:,idx(i)),handles.Istack(:,:,template_idx),transform,optimizer,metric);
    waitbar(i/length(idx),hwait,'Please wait while aligning image frames to template');
end
handles.Istack_warped = Istack_warped;
msgbox('Image registration complete. Use "Pre-processing -> Motion correction -> View corrected stack" to ensure motion artifact is removed. Once satisfied, save the updated image stack');
try
    delete(hwait);
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function save_stack_Callback(hObject, eventdata, handles)
% hObject    handle to save_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
if(~isfield(handles,'Istack_warped') || isempty(handles.Istack_warped))
    
    errordlg('You must first align multiple image frames to a template','Bad input','modal');
    return;
end

answer = questdlg(sprintf('You are about to overwrite the original image stack with motion-corrected stack. Are you sure you want to proceed? Click "No" to go back and view the results of image registration'),'Proceed?','Yes','No','No');
switch answer
    case 'No'
        return;
    case 'Yes'
        hwait = waitbar(0,sprintf('Please wait while overwriting %s with motion-corrected data',handles.curr_file));
        imwrite(handles.Istack_warped(:,:,1),fullfile(handles.curr_folder,handles.curr_file));
        for i=2:size(handles.Istack_warped,3)
            if(~mod(i,100))
                waitbar(i/size(handles.Istack_warped,3),hwait,sprintf('Please wait while overwriting %s with motion-corrected data',handles.curr_file));
            end
            imwrite(handles.Istack_warped(:,:,i),fullfile(handles.curr_folder,handles.curr_file),'writemode','append');
        end
        handles.Istack = handles.Istack_warped;
end
try
    delete(hwait);
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function view_corrected_stack_Callback(hObject, eventdata, handles)
% hObject    handle to view_corrected_stack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
if(~isfield(handles,'Istack_warped') || isempty(handles.Istack_warped))
    
    errordlg('You must first align multiple image frames to a template','Bad input','modal');
    return;
end
ViewImageStack(handles.Istack_warped);


% --------------------------------------------------------------------
function correlation_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to correlation_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    handles.curr_folder_idx;
    handles.curr_folder;
    handles.curr_file;
catch
    errordlg('You must first select a file','Bad input','modal');
    return;
end
load('params.mat');
hwait = waitbar(0,'Please wait while computing pair-wise image frame correlations.',...
    'CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');
setappdata(hwait,'canceling',0)
% Read the image stack
if(isfield(handles,'Istack') && ~isempty(handles.Istack))
    ;
else
    multiWaitbar('CloseAll','Name','');
    multiWaitbar('Reading image stack.','Busy','Name','Please wait');
    handles.Istack = ReadTiffStack(fullfile(handles.curr_folder,handles.curr_file));
    multiWaitbar('CloseAll','Name','');
end

% Compute pair-wise correlations
N = size(handles.Istack,3);
C = zeros(N,N);

for i=1:N
    if getappdata(hwait,'canceling')
        break
    end
    waitbar(i/N,hwait,'Please wait while computing pair-wise image frame correlations');
    if(params.parallel)
        parfor j=i:N
            C(i,j) = corr2(handles.Istack(:,:,i),handles.Istack(:,:,j));
        end
    else
        for j=i:N
            C(i,j) = corr2(handles.Istack(:,:,i),handles.Istack(:,:,j));
        end
    end
end
for i=1:N
    for j=i:N
        C(j,i) = C(i,j);
    end
end
try
    delete(hwait);
end
figure;
imagesc(C); colorbar; axis square;
colormap('jet');
title(sprintf('Pair-wise image frame correlations. Low values indicate possible misalignment (motion artifact)'));
xlabel('Frame #'); ylabel('Frame #');


% --------------------------------------------------------------------
function dF_Callback(hObject, eventdata, handles)
% hObject    handle to dF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
figure;
t = 0:1/data.fps:size(data.F_cell,2)/data.fps-1/data.fps;
sig = data.dF_cell;

% calculate shift
mi = min(sig,[],2);
ma = max(sig,[],2);
shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,1,size(sig,2));

%plot data
plot(t,sig+shift)

% edit axes
set(gca,'ytick',mean(sig+shift,2),'yticklabel',1:data.N)
grid on
ylim([mi(1) max(max(shift+sig))]);
xlabel('Time (s)'); ylabel('Neuron ID');
title(data.filename);

% --------------------------------------------------------------------
function recompute_Callback(hObject, eventdata, handles)
% hObject    handle to recompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    errordlg('You must first select either "Process this file" or "Batch process" then use this tool to recompute synchronization and functional connectivity if you made changes to the event detection or preferences','Bad input','modal');
    return;
end
try
    load('params.mat');
catch
    errordlg('Could not load params.mat. Please run Analyze -> Preferences -> Revert to defaults','File not found','modal');
    return;
end
handles.data.params = params;
[~,y]=fileparts(handles.data.filename);
multiWaitbar('CloseAll','Name','');
multiWaitbar('Synchronization Cluster Analysis',0,'Name',['Recomputing analysis for ' y]);
if(params.FC.method_idx==1)
    multiWaitbar('Functional connectivity: cross-correlation',0,'Name',y);
elseif(params.FC.method_idx==2)
    multiWaitbar('Functional connectivity: partial correlation',0,'Name',y);
elseif(params.FC.method_idx==3)
    multiWaitbar('Functional connectivity: phase',0,'Name',y);
elseif(params.FC.method_idx==4)
    multiWaitbar('Functional connectivity: Granger causality',0,'Name',y);
elseif(params.FC.method_idx==5)
    multiWaitbar('Functional connectivity: Transfer entropy',0,'Name',y);
elseif(params.FC.method_idx==6)
    multiWaitbar('Functional connectivity: cross-correlation',0,'Name',y);
    multiWaitbar('Functional connectivity: partial correlation',0,'Name',y);
    multiWaitbar('Functional connectivity: phase',0,'Name',y);
    multiWaitbar('Functional connectivity: Granger causality',0,'Name',y);
    multiWaitbar('Functional connectivity: Transfer entropy',0,'Name',y);
end
multiWaitbar('Detecting network ensembles',0);
multiWaitbar('Computing calcium event transient kinetics',0);
multiWaitbar('Inferring spikes from fluorescence',0);

multiWaitbar('Saving results and making summary figure',0);
C = SCA(handles.data,'DISPLAY',0,'save_flag',0,'wb',0);
multiWaitbar('Synchronization Cluster Analysis',1,'Color','g');
SI = max(C.SI);
handles.data.SynchroCluster = C;
handles.data.phase = C.phase;
handles.data.SI = SI;
N = handles.data.N;
if(params.FC.method_idx==1)
    [A,C] = FC_crosscorr(handles.data);
    clu = clustering_coef_bu(A);
    
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.CC.A = A;
    handles.data.FC.CC.C = C;
    handles.data.FC.CC.clustering_coef = clu;
    handles.data.FC.CC.modularity_Ci = Ci;
    handles.data.FC.CC.modularity_Q = Q;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.CC.Controllability.Nd = Nd;
    handles.data.FC.CC.Controllability.drivernodes = drivernodes;
    handles.data.FC.CC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: cross-correlation',1,'Color','g');
    
elseif(params.FC.method_idx==2)
    multiWaitbar('Functional connectivity: partial correlation','Busy');
    [A,rho] = FC_partialcorr(handles.data.dF_cell);
    clu = clustering_coef_bu(A);
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.PC.A = A;
    handles.data.FC.PC.rho = rho;
    handles.data.FC.PC.clustering_coef = clu;
    handles.data.FC.PC.modularity_Ci = Ci;
    handles.data.FC.PC.modularity_Q = Q;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.PC.Controllability.Nd = Nd;
    handles.data.FC.PC.Controllability.drivernodes = drivernodes;
    handles.data.FC.PC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: partial correlation',1,'Color','g');
elseif(params.FC.method_idx==3)
    [A,P] = FC_phase(handles.data);
    clu = clustering_coef_bu(A);
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.phase.A=A;
    handles.data.FC.phase.Pvals = P;
    handles.data.FC.phase.clustering_coef = clu;
    handles.data.FC.phase.modularity_Ci = Ci;
    handles.data.FC.phase.modularity_Q = Q;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.phase.Controllability.Nd = Nd;
    handles.data.FC.phase.Controllability.drivernodes = drivernodes;
    handles.data.FC.phase.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: phase',1,'Color','g');
elseif(params.FC.method_idx==4)
    [A,P,F] = FC_granger(handles.data.dF_cell);
    clu = clustering_coef_bd(A);
    try
        [Ci,Q] = modularity_louvain_dir(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    [~,~,d] = degrees_dir(A);
    handles.data.FC.GC.A = A;
    handles.data.FC.GC.Pvals = P;
    handles.data.FC.GC.F = F;
    handles.data.FC.GC.clustering_coef = clu;
    handles.data.FC.GC.modularity_Ci = Ci;
    handles.data.FC.GC.modularity_Q = Q;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.GC.Controllability.Nd = Nd;
    handles.data.FC.GC.Controllability.drivernodes = drivernodes;
    handles.data.FC.GC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: Granger causality',1,'Color','g');
elseif(params.FC.method_idx==5)
    [peakTE,CI] = FC_transfer_entropy(handles.data);
    A = peakTE;
    clu = clustering_coef_wd(peakTE);
    try
        [Ci,Q] = modularity_louvain_dir(peakTE);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    [~,~,d] = degrees_dir(peakTE);
    handles.data.FC.TE.peakTE = peakTE;
    handles.data.FC.TE.A = peakTE;
    handles.data.FC.TE.CI = CI;
    handles.data.FC.TE.clustering_coef = clu;
    handles.data.FC.TE.modularity_Ci = Ci;
    handles.data.FC.TE.modularity_Q = Q;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.TE.Controllability.Nd = Nd;
    handles.data.FC.TE.Controllability.drivernodes = drivernodes;
    handles.data.FC.TE.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: Transfer entropy',1,'Color','g');
elseif(params.FC.method_idx==6)
    % Do all
    [A,C] = FC_crosscorr(handles.data);
    clu = clustering_coef_bu(A);
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.CC.modularity_Ci = Ci;
    handles.data.FC.CC.modularity_Q = Q;
    handles.data.FC.CC.A = A;
    handles.data.FC.CC.C = C;
    handles.data.FC.CC.clustering_coef = clu;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.CC.Controllability.Nd = Nd;
    handles.data.FC.CC.Controllability.drivernodes = drivernodes;
    handles.data.FC.CC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: cross-correlation',1,'Color','g');
    
    multiWaitbar('Functional connectivity: partial correlation','Busy');
    [A,rho] = FC_partialcorr(handles.data.dF_cell);
    clu = clustering_coef_bu(A);
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.PC.modularity_Ci = Ci;
    handles.data.FC.PC.modularity_Q = Q;
    handles.data.FC.PC.clustering_coef = clu;
    handles.data.FC.PC.A = A;
    handles.data.FC.PC.rho = rho;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.PC.Controllability.Nd = Nd;
    handles.data.FC.PC.Controllability.drivernodes = drivernodes;
    handles.data.FC.PC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: partial correlation',1,'Color','g');
    
    [A,P] = FC_phase(handles.data);
    clu = clustering_coef_bu(A);
    try
        [Ci,Q] = modularity_louvain_und(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    d = sum(A);
    handles.data.FC.phase.modularity_Ci = Ci;
    handles.data.FC.phase.modularity_Q = Q;
    handles.data.FC.phase.A=A;
    handles.data.FC.phase.Pvals = P;
    handles.data.FC.phase.clustering_coef = clu;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.phase.Controllability.Nd = Nd;
    handles.data.FC.phase.Controllability.drivernodes = drivernodes;
    handles.data.FC.phase.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: phase',1,'Color','g');
    
    [A,P,F] = FC_granger(handles.data.dF_cell);
    clu = clustering_coef_bd(A);
    try
        [Ci,Q] = modularity_louvain_dir(A);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    [~,~,d] = degrees_dir(A);
    handles.data.FC.GC.modularity_Ci = Ci;
    handles.data.FC.GC.modularity_Q = Q;
    handles.data.FC.GC.clustering_coef = clu;
    handles.data.FC.GC.A = A;
    handles.data.FC.GC.Pvals = P;
    handles.data.FC.GC.F = F;
    [Nd,drivernodes,Nconfigs] = ExactControllability(A,'plotting',0);
    handles.data.FC.GC.Controllability.Nd = Nd;
    handles.data.FC.GC.Controllability.drivernodes = drivernodes;
    handles.data.FC.GC.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: Granger causality',1,'Color','g');
    
    [peakTE,CI] = FC_transfer_entropy(handles.data);
    clu = clustering_coef_wd(peakTE);
    try
        [Ci,Q] = modularity_louvain_dir(peakTE);
    catch
        Ci = zeros(N,1);
        Q = 0;
    end
    [~,~,d] = degrees_dir(peakTE);
    handles.data.FC.TE.modularity_Ci = Ci;
    handles.data.FC.TE.modularity_Q = Q;
    handles.data.FC.TE.peakTE = peakTE;
    handles.data.FC.TE.A = peakTE;
    handles.data.FC.TE.CI = CI;
    handles.data.FC.TE.clustering_coef = clu;
    [Nd,drivernodes,Nconfigs] = ExactControllability(peakTE,'plotting',0);
    handles.data.FC.TE.Controllability.Nd = Nd;
    handles.data.FC.TE.Controllability.drivernodes = drivernodes;
    handles.data.FC.TE.Controllability.Nconfigs = Nconfigs;
    multiWaitbar('Functional connectivity: Transfer entropy',1,'Color','g');
    
end
if(params.FC.method_idx==6) % Default to phase method if the user selects "all" for the purpose of making figure
    A = handles.data.FC.phase.A;
    clu = handles.data.FC.phase.clustering_coef;
end

% Sychronization by modules
SI_m = zeros(max(Ci),1);
for k=1:max(Ci)
    tmp_s.F_cell = handles.data.F_cell(Ci==k,:);
    tmp_s.Spikes_cell = handles.data.Spikes_cell(Ci==k);
    tmp_s.fps = handles.data.fps;
    c_mod = SCA(tmp_s);
    SI_m(k) = max(c_mod.SI);
end

% Mean oscillation period (s)

ISI = cell(handles.data.N,1);
for k=1:N
    ISI{k} = diff(handles.data.Spikes_cell{k});
end
OP = zeros(N,1);
for k=1:N
    try
        OP(k) = mean(ISI{k})/handles.data.fps;
    end
end
% Network ensembles
[Nhigh_activity_frames,ensemble_per_second,ensemble_frames,CoreEnsembles,EnsembleR]=NetworkEnsemble(handles.data);
multiWaitbar('Detecting network ensembles',1,'Color','g');
handles.data.NetworkEnsemble.ensembles_per_second = ensemble_per_second;
handles.data.NetworkEnsemble.ensemble_frames = ensemble_frames;
handles.data.NetworkEnsemble.Nensembles = Nhigh_activity_frames;
handles.data.NetworkEnsemble.CoreEnsembles = CoreEnsembles;
handles.data.NetworkEnsemble.CorrelatedEnsembles = EnsembleR;
fprintf('\tDetermining transient kinetics - amplitude, rise time, and fall time\n');

[DF,rise_time,fall_time,CV] = GetTransients(handles.data);
multiWaitbar('Computing calcium event transient kinetics',1,'Color','g');
ci = bootci(40,@(x) mean(x),DF(~isnan(DF)));

Rise = [rise_time{:}]; Fall = [fall_time{:}];
ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));

ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
fprintf('\tPerforming spike inference\n');
% Fast spike inference
spk = zeros(size(handles.data.F_cell));

for m=1:handles.data.N
    multiWaitbar('Inferring spikes from fluorescence',m/handles.data.N);
    x = run_oopsi(handles.data.dF_cell(m,:));
    spk(m,:) = x.n';
    spk(m,:) = spk(m,:)./max(spk(m,:));
end
multiWaitbar('Inferring spikes from fluorescence',1,'Color','g');
multiWaitbar('Saving results and making summary figure','Busy');


handles.data.SI_m = SI_m;
handles.data.DF = DF;
handles.data.rise_time = rise_time;
handles.data.fall_time = fall_time;
handles.data.CV = CV;
handles.data.modules = Ci;
handles.data.modularity = Q;
handles.data.foopsi = spk;
fid = fopen([handles.data.filename(1:end-4) '_summary.txt'],'w');
data = handles.data;

C = data.SynchroCluster;
N = data.N;

dat = zeros(N,10);
OP = zeros(N,1);
for k=1:N
    try
        OP(k) = mean(ISI{k})/data.fps;
    end
end
for k=1:N
    x = sort(data.F_cell(k,:));
    dat(k,1) = mean(x(1:ceil(.1*length(x))));
    dat(k,2) = numel(data.Spikes_cell{k});
    dat(k,3) = OP(k);
    dat(k,4) = data.DF(k);
    dat(k,5) = data.CV(k);
    r = data.rise_time{k}; tau = data.fall_time{k};
    dat(k,6) = mean(r(~isnan(r)));
    dat(k,7) = mean(tau(~isnan(tau)));
    dat(k,8) = sum(A(k,:))/N;
    dat(k,9) = max(C.PI(k,:));
    dat(k,10) = nnz(C.PI(k,:)>.01);
    dat(k,11) = clu(k);
end

% Export everything to .txt file
% Mean oscillation period (s)
ISI = cell(N,1);
for k=1:N
    ISI{k} = diff(data.Spikes_cell{k});
end


ci = bootci(40,@(x) mean(x),data.DF(~isnan(data.DF)));

Rise = [data.rise_time{:}]; Fall = [data.fall_time{:}];
ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));

ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
if(fid)
    fprintf(fid,'Summary for file:\t %s\n\n',data.filename);
    fprintf(fid,'Total cells (ROIs):\t %d\n',N);
    fprintf(fid,'Frames:\t %d\nframe rate:\t %d\nduration (s):\t %.02f\n',...
        size(data.F_cell,2),data.fps,...
        size(data.F_cell,2)/data.fps);
    % Output functional connectivity summaries
    if(params.FC.method_idx==1)
        fprintf(fid,'Functional connectivity method: cross-correlation\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.CC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.CC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.CC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.CC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.CC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.CC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.CC.Controllability.Nconfigs);
    end
    
    if(params.FC.method_idx==2)
        fprintf(fid,'Functional connectivity method: partial-correlation\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.PC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.PC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.PC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.PC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.PC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.PC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.PC.Controllability.Nconfigs);
    end
    if(params.FC.method_idx==3)
        fprintf(fid,'Functional connectivity method: phase\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.phase.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.phase.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.phase.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.phase.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.phase.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.phase.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.phase.Controllability.Nconfigs);
    end
    
    if(params.FC.method_idx==4)
        fprintf(fid,'Functional connectivity method: Granger causality\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.GC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.GC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.GC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.GC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.GC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.GC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.GC.Controllability.Nconfigs);
    end
    if(params.FC.method_idx==5)
        fprintf(fid,'Functional connectivity method: Transfer entropy\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.TE.peakTE))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.TE.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.TE.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.TE.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.TE.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.TE.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.TE.Controllability.Nconfigs);
    end
    if(params.FC.method_idx==6)
        fprintf(fid,'Functional connectivity method: cross-correlation\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.CC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.CC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.CC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.CC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.CC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.CC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.CC.Controllability.Nconfigs);
        
        fprintf(fid,'Functional connectivity method: partial-correlation\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.PC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.PC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.PC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.PC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.PC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.PC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.PC.Controllability.Nconfigs);
        fprintf(fid,'Functional connectivity method: phase\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.phase.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.phase.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.phase.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.phase.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.phase.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.phase.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.phase.Controllability.Nconfigs);
        fprintf(fid,'Functional connectivity method: Granger causality\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.GC.A))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.GC.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.GC.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.GC.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.GC.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.GC.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.GC.Controllability.Nconfigs);
        fprintf(fid,'Functional connectivity method: Transfer entropy\n');
        fprintf(fid,'\tMean global connectivity:\t %.04f\n',mean(sum(data.FC.TE.peakTE))/N);
        
        fprintf(fid,'\tModularity:\t %.02f\n',data.FC.TE.modularity_Q);
        fprintf(fid,'\tNumber of modules:\t %d\n',max(data.FC.TE.modularity_Ci));
        fprintf(fid,'\tNumber of driver nodes:\t %d\n',data.FC.TE.Controllability.Nd);
        fprintf(fid,'\tList of driver nodes:\t');
        for q=1:length(data.FC.TE.Controllability.drivernodes)
            fprintf(fid,'%d, ',data.FC.TE.Controllability.drivernodes(q));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\tNumber of possible driver node configurations:\t %d\n',data.FC.TE.Controllability.Nconfigs);
        
    end
    
    fprintf(fid,'Global synchronization index:\t %.02f\n',SI);
    fprintf(fid,'Synchronization index by modules:\t ');
    for k=1:length(SI_m)
        fprintf(fid,'%f\t',SI_m(k));
    end
    fprintf(fid,'\n');
    
    fprintf(fid,'Number of network ensembles:\t %d\n',data.NetworkEnsemble.Nensembles);
    fprintf(fid,'Ensembles per second:\t %.04f\n',data.NetworkEnsemble.ensembles_per_second);
    for q=1:length(data.NetworkEnsemble.CoreEnsembles)
        fprintf(fid,'Core ensemble %d:\t',q);
        ens = data.NetworkEnsemble.CoreEnsembles{q};
        for w=1:length(ens)
            fprintf(fid,'%d, ',ens(w));
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'Mean oscillation period (s):\t %.02f\n',mean(OP));
    fprintf(fid,'Mean amplitude (deltaF/F0):\t %.04f\n95%% CI:\t%.04f\t %.04f\n',mean(DF(~isnan(DF))),ci(1),ci(2));
    fprintf(fid,'Coefficient of variation in amplitude:\t%f\n',mean(CV(~isnan(CV))));
    fprintf(fid,'Mean rise time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2));
    fprintf(fid,'Mean fall time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2));
    %Now output cell based report
    fprintf(fid,'\n\nNeuron ID\tBaseline fluorescence\tTotal events\t<ISI> (s)\tAmplitude\tCV\t<Rise time> (s)\t<Fall time> (s)\tConnectivity\tParticipation\t# assemblies\tClustering coef\n');
    for k=1:data.N
        fprintf(fid,'%d\t%f\t%.02f\t%.04f\t%.02f\t%f\t%.02f\t%.04f\t%.04f\t%d\t%.04f\t%0.04f\n',k,dat(k,1),dat(k,2),dat(k,3),dat(k,4),dat(k,5),dat(k,6),dat(k,7),dat(k,8),dat(k,9),dat(k,10),dat(k,11));
    end
    fclose(fid);
    
    cprintf('*blue','%s\n', ['Summary statistics written to ' data.filename(1:end-4) '_summary.txt']);
else
    cprintf('*red','%s\n','ERROR: Could not write summary statistics. Unable to create file.');
    
end

% Finally, make a composite figure and save it as eps
if(ispc)
    slash = '\';
else
    slash = '/';
end
netfig = figure;
set (netfig, 'Units', 'normalized', 'Position', [0,0,1,1],'PaperOrientation','portrait','PaperType','arch-C');
subplot(2,2,1);
if(isfield(handles.data,'image'))
    if(size(handles.data.image,3)~=3)
        image = zeros([size(handles.data.image) 3]);
        image(:,:,2) = imadjust(handles.data.image);
    elseif(size(handles.data.image,3)==3)
        image = handles.data.image;
    else
        image = logical(handles.data.L);
    end
end
imagesc(image);
hold on

c = regionprops(handles.data.L,'Centroid');
c = reshape([c.Centroid],2,handles.data.N);
x = c(1,:);
y = c(2,:);
Colors = varycolor(max(Ci));
if(isempty(Colors))
    Colors = [0 0 1];
end
gplot(A,[x' y']);
set(findobj(gcf,'Type','Line'),'LineStyle',':','LineWidth',.01,'MarkerSize',6,'Marker','o','MarkerFaceColor',[1 1 0],'Color',[1 0 0]);

plot(x(1,~logical(sum(A)) & ~logical(sum(A'))),y(1,~logical(sum(A)) & ~logical(sum(A'))),'ro','MarkerFaceColor',[1 0 0])
set(gca,'YDir','reverse');
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
title({['N = ' num2str(handles.data.N) ',Dropped = ' num2str(sum(~logical(sum(A)) & ~logical(sum(A'))))]});

subplot(2,2,2);
[~,pos] = max(C.PI,[],2);
[~,IDX] = sort(pos);
imagesc(C.C(IDX,IDX),[0 1]); axis square; colorbar;
colormap('jet');
xlabel('Neuron ID (rearranged)'); ylabel('Neuron ID (rearranged)');
msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
title(msg);

subplot(2,2,3);

set(gca, 'ColorOrder', Colors);
hold all;
if(sum(Ci)==0)
    Ci = ones(N,1);
    scatter(x,y,'filled','SizeData',10^2);
else
    for i=1:max(Ci)
        scatter(x(Ci==i),y(Ci==i),'filled','SizeData',10^2);
    end
end
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(gca,'YDir','reverse');

title_string = sprintf('N = %d. Modules = %d. Modularity = %f\n',length(x),max(Ci),Q);
title_string = sprintf('%sSynchronization by modules: ',title_string);
for k=1:length(SI_m)
    title_string = sprintf('%s %f,',title_string,SI_m(k));
end
title(title_string);

subplot(2,2,4);
t = 0:1/handles.data.fps:size(handles.data.F_cell,2)/handles.data.fps-1/handles.data.fps;
hold all;
cntr = 1;
for q=1:max(Ci)
    modules = find(Ci==q);
    for m=1:length(modules)
        try
            evnts = handles.data.Spikes_cell{modules(m)};
            plot(evnts/handles.data.fps,cntr,'.','Color',Colors(q,:));
            cntr = cntr+1;
        end
    end
end
title('Activity rearranged by modules');
frames = size(handles.data.F_cell,2);
ylim([0 N]); xlim([0 frames/handles.data.fps]);
xlabel('Time (s)'); ylabel('Neuron ID (rearranged)');
% Save this figure in a new folder
[pathstr, name] = fileparts(handles.data.filename);
if(~exist([pathstr slash 'Figures'],'dir'))
    mkdir(pathstr,'Figures');
end
set(netfig,'PaperPositionMode','auto');
print(netfig,'-depsc',[pathstr slash 'Figures' slash name '.eps']);

print(netfig,'-dtiff',[pathstr slash 'Figures' slash name '.tif']);
close(netfig);
drawnow;


save([pathstr slash 'analysis-' name '.mat'],'data');
multiWaitbar('Saving results and making summary figure',1,'Color','g');
pause(1);
multiWaitbar('CloseAll','Name','');


% --------------------------------------------------------------------
function synchronization_matrix_Callback(hObject, eventdata, handles)
% hObject    handle to synchronization_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
handles.data = data;
C = handles.data.SynchroCluster;
figure;
imagesc(C.C); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID'); ylabel('Neuron ID');
msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
title(msg);
figure;

[~,pos] = max(C.PI,[],2);
[~,IDX] = sort(pos);
imagesc(C.C(IDX,IDX),[0 1]); axis square; colorbar;
colormap('jet');
xlabel('Neuron ID (rearranged)'); ylabel('Neuron ID (rearranged)');
msg = sprintf('Global synchronization index: %f\n%d clusters detected',max(C.SI),size(C.PI,2));
title(msg);


% --------------------------------------------------------------------
function Untitled_5_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function export_raw_fluorescence_Callback(hObject, eventdata, handles)
% hObject    handle to export_raw_fluorescence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save raw fluorescence data as:');
dlmwrite(fullfile(pathname,filename),data.F_cell);
icon = imread('480.png');
msgbox('Raw fluorescence data saved to txt file. Each row is an ROI and each column is a frame; value at row i, column j is the mean fluorescence intensity of ROI i in frame j.','Saved','custom',icon);


% --------------------------------------------------------------------
function export_dF_Callback(hObject, eventdata, handles)
% hObject    handle to export_dF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save deltaF/F data as:');
dlmwrite(fullfile(pathname,filename),data.dF_cell);
icon = imread('480.png');
msgbox('deltaF/F saved to txt file. Each row is an ROI and each column is a frame; value at row i, column j is the deltaF/F of ROI i in frame j.','Saved','custom',icon);

% --------------------------------------------------------------------
function export_spikes_Callback(hObject, eventdata, handles)
% hObject    handle to export_spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save inferred spikes as:');
dlmwrite(fullfile(pathname,filename),data.foopsi);
icon = imread('480.png');
msgbox('Inferred spikes saved to txt file. Each row is an ROI and each column is a frame; value at row i, column j is the normalized spike probability of ROI i in frame j.','Saved','custom',icon);

% --------------------------------------------------------------------
function export_transients_Callback(hObject, eventdata, handles)
% hObject    handle to export_transients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save onsets of calcium transients as:');
try
    fid = fopen(fullfile(pathname,filename),'w');
catch
    errordlg(['Could not open ' fullfile(pathname,filename) ' for writting'],'modal');
    return;
end
for i=1:data.N
    fprintf(fid,'%d: ',i);
    evts = data.Spikes_cell{i};
    for j=1:length(evts)
        fprintf(fid,'%d,',evts(j));
    end
    fprintf(fid,'\n');
end
fclose(fid);


icon = imread('480.png');
msgbox('Onsets of calcium transients saved to txt file. File format is "ROI #: frame1,frame2,..."','Saved','custom',icon);


% --------------------------------------------------------------------
function export_synchronization_Callback(hObject, eventdata, handles)
% hObject    handle to export_synchronization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save synchronization matrix as:');
dlmwrite(fullfile(pathname,filename),data.SynchroCluster.C);
icon = imread('480.png');
msgbox('Pair-wise synchronization matrix saved to txt file.','Saved','custom',icon);


% --------------------------------------------------------------------
function export_FC_Callback(hObject, eventdata, handles)
% hObject    handle to export_FC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save functional connectivity adjacency matrix as:');
dlmwrite(fullfile(pathname,filename),data.SynchroCluster.A);
icon = imread('480.png');
msgbox('Functional connectivity adjacency matrix saved to txt file.','Saved','custom',icon);


% --------------------------------------------------------------------
function export_kinetics_Callback(hObject, eventdata, handles)
% hObject    handle to export_kinetics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
        errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
        return;
    end
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
catch
    errordlg(['Please select a file first and make sure an analysis-' handles.curr_file(1:end-4) '.mat exists'],'modal');
    return;
end
[filename,pathname] = uiputfile('.txt','Save single-cell summaries as:');
try
    fid = fopen(fullfile(pathname,filename),'w');
catch
    errordlg(['Could not open ' fullfile(pathname,filename) ' for writting'],'modal');
    return;
end
C = data.SynchroCluster;
N = data.N;

dat = zeros(N,10);
clu = clustering_coef_bu(C.A);
for k=1:N
    x = sort(data.F_cell(k,:));
    dat(k,1) = mean(x(1:ceil(.1*length(x))));
    dat(k,2) = numel(data.Spikes_cell{k});
    dat(k,3) = mean(diff(data.Spikes_cell{k}))/data.fps;
    dat(k,4) = data.DF(k);
    dat(k,5) = data.CV(k);
    r = data.rise_time{k}; tau = data.fall_time{k};
    dat(k,6) = mean(r(~isnan(r)));
    dat(k,7) = mean(tau(~isnan(tau)));
    dat(k,8) = sum(C.A(k,:))/N;
    dat(k,9) = max(C.PI(k,:));
    dat(k,10) = nnz(C.PI(k,:)>.01);
    dat(k,11) = clu(k);
end

% Export everything to .txt file
% Mean oscillation period (s)
ISI = cell(N,1);
for k=1:N
    ISI{k} = diff(data.Spikes_cell{k});
end

OP = mean([ISI{:}])/data.fps;
ci = bootci(40,@(x) mean(x),data.DF(~isnan(data.DF)));

Rise = [data.rise_time{:}]; Fall = [data.fall_time{:}];
ci_r = bootci(40,@(x) mean(x),Rise(~isnan(Rise)));

ci_tau = bootci(40,@(x) mean(x),Fall(~isnan(Fall)));
if(fid)
    fprintf(fid,'Summary for file:\t %s\n\n',data.filename);
    fprintf(fid,'Total cells (ROIs):\t %d\n',N);
    fprintf(fid,'Frames:\t %d\nframe rate:\t %d\nduration (s):\t %.02f\n',...
        size(data.F_cell,2),data.fps,...
        size(data.F_cell,2)/data.fps);
    fprintf(fid,'Modularity:\t %.04f\nNumber of modules:\t %d\n',data.modularity,max(data.modules(:)));
    fprintf(fid,'Global synchronization index:\t %.02f\n',data.SI);
    fprintf(fid,'Synchronization index by modules:\t ');
    for k=1:length(data.SI_m)
        fprintf(fid,'%f\t',data.SI_m(k));
    end
    fprintf(fid,'\n');
    fprintf(fid,'Mean global functional connectivity:\t %.04f\n',mean(sum(C.A))/data.N);
    fprintf(fid,'Mean oscillation period (s):\t %.02f\n',OP);
    fprintf(fid,'Mean amplitude (deltaF/F0):\t %.04f\n95%% CI:\t%.04f\t %.04f\n',mean(data.DF(~isnan(data.DF))),ci(1),ci(2));
    fprintf(fid,'Coefficient of variation in amplitude:\t%f\n',mean(data.CV(~isnan(data.CV))));
    fprintf(fid,'Mean rise time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Rise(~isnan(Rise))),ci_r(1),ci_r(2));
    fprintf(fid,'Mean fall time (s):\t %.02f\n95%% CI:\t%.02f\t%.02f\n',mean(Fall(~isnan(Fall))),ci_tau(1),ci_tau(2));
    %Now output cell based report
    fprintf(fid,'\n\nNeuron ID\tBaseline fluorescence\tTotal events\t<ISI> (s)\tAmplitude\tCV\t<Rise time> (s)\t<Fall time> (s)\tConnectivity\tParticipation\t# assemblies\tClustering coef\n');
    for k=1:data.N
        fprintf(fid,'%d\t%f\t%.02f\t%.04f\t%.02f\t%f\t%.02f\t%.04f\t%.04f\t%d\t%.04f\t%0.04f\n',k,dat(k,1),dat(k,2),dat(k,3),dat(k,4),dat(k,5),dat(k,6),dat(k,7),dat(k,8),dat(k,9),dat(k,10),dat(k,11));
    end
    fclose(fid);
end
icon = imread('480.png');
msgbox('Single cell summaries saved to txt file.','Saved','custom',icon);


% --------------------------------------------------------------------
function preferences_Callback(hObject, eventdata, handles)
% hObject    handle to preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg = sprintf('Select transform type for image registration.\n1 = translation, 2 = rigid, 3 = similarity, 4 = affine');
try
    answer = inputdlg({msg,'Number of iterations'},'Image registration preferences',1,{'1','100'});
    
    handles.transformtype = str2num(answer{1});
    handles.transforiter = str2num(answer{2});
catch
    errordlg('Bad input format. Please try again','Bad input','modal');
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FC_CC_Callback(hObject, eventdata, handles)
% hObject    handle to FC_CC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.FC.CC;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used cross-correlation as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end

% Show cross-correlation matrix and binary adjacency matrix
figure;
subplot(1,2,1); imagesc(handles.data.FC.CC.C); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Pairwise cross-correlation matrix.\nEntries are max normalized Pearson coefficient between [-lag,lag]'));
subplot(1,2,2); imagesc(handles.data.FC.CC.A); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Binary adjacency matrix\nN=%d surrogate resampling and max lag %.2f (s)',handles.data.params.FC.CC.Nsur,handles.data.params.FC.CC.maxlag));

% --------------------------------------------------------------------
function FC_PC_Callback(hObject, eventdata, handles)
% hObject    handle to FC_PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.FC.PC;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used partial correlation as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end

% Show partial correlation matrix and binary adjacency matrix
figure;
subplot(1,2,1); imagesc(handles.data.FC.PC.rho); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Pairwise partial correlation matrix.\nEntries are normalized partial correlation coefficients'));
subplot(1,2,2); imagesc(handles.data.FC.PC.A); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Binary adjacency matrix\nThresholding at alpha level %f for significance',handles.data.params.FC.PC.alpha));

% --------------------------------------------------------------------
function FC_phase_Callback(hObject, eventdata, handles)
% hObject    handle to FC_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.FC.phase;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used "phase" as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end

% Show pair-wise phase similarity matrix and binary adjacency matrix
figure;
subplot(1,2,1); imagesc(handles.data.SynchroCluster.C); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Pairwise phase-locking matrix\nEntries are circular variance of pair-wise phase difference'));
subplot(1,2,2); imagesc(handles.data.FC.phase.A); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Binary adjacency matrix\nN=%d resampling, alpha level %f',handles.data.params.FC.phase.Nsur,handles.data.params.FC.phase.alpha));

% --------------------------------------------------------------------
function FC_GC_Callback(hObject, eventdata, handles)
% hObject    handle to FC_GC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.FC.GC;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used Granger causality as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
figure;
subplot(1,2,1); imagesc(handles.data.FC.GC.F); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title('Pairwise Granger causality F statistic');
subplot(1,2,2); imagesc(handles.data.FC.GC.A); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Binary adjacency matrix\nThresholding at alpha level %f for significance',handles.data.params.FC.GC.alpha));

% --------------------------------------------------------------------
function FC_TE_Callback(hObject, eventdata, handles)
% hObject    handle to FC_TE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.FC.TE;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
figure;
subplot(1,2,1); imagesc(handles.data.FC.TE.peakTE); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Pairwise max transfer entropy from lags %d to %d',handles.data.params.FC.TE.lags(1),handles.data.params.FC.TE.lags(end)));
subplot(1,2,2); imagesc(handles.data.FC.TE.CI); colorbar; axis square;
colormap('jet');
xlabel('Neuron ID (to)'); ylabel('Neuron ID (from)');
title(sprintf('Pairwise coincidence index from lags %d to %d',handles.data.params.FC.TE.lags(1),handles.data.params.FC.TE.lags(end)));


% --------------------------------------------------------------------
function network_ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to network_ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg1 = sprintf('Number of standard deviations above 0 to binarize spike probabilities');
msg2 = sprintf('Number of times to reshuffle binary activity data for statistical significance:');

try
    load params.mat
end
try
    answer = inputdlg({msg1,msg2},'Network ensemble',1,...
        {num2str(params.network_ensemble_sd),num2str(params.network_ensemble_Nsur)});
catch
    answer = inputdlg({msg1,msg2},'Network ensemble',1,{'3','1000'});
end
params.network_ensemble_sd = str2num(answer{1});
params.network_ensemble_Nsur = str2num(answer{2});

save('FluoroSNNAP_code/params.mat','params');


% --------------------------------------------------------------------
function ensemble_Callback(hObject, eventdata, handles)
% hObject    handle to ensemble (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    handles.data.NetworkEnsemble;
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
spk = data.foopsi;

% Binarize spks. Prompt user for number of s.d above 0

APs = zeros(size(spk));
for i=1:size(spk,1)
    sd = std(spk(i,:));
    APs(i,:) = spk(i,:)>sd*data.params.network_ensemble_sd;
end
simult = sum(APs,1);
t = 0:1/data.fps:size(data.F_cell,2)/data.fps-1/data.fps;
figure
plot(t,simult./data.N*100,'k','LineWidth',2); xlabel('Time (s)'); ylabel('% of cells co-active');
figure
hold on
for i=1:size(APs,2)
    idx = find(APs(:,i));
    if(~isempty(idx))
        if(nnz(find(i==data.NetworkEnsemble.ensemble_frames)))
            plot(t(i),idx,'r.');
        else
            plot(t(i),idx,'b.');
        end
    end
end
xlabel('Time (s)'); ylabel('Neuron ID');
title(sprintf('Network ensemles in red\n%.04f ensembles per second, %d total ensembles',data.NetworkEnsemble.ensembles_per_second,data.NetworkEnsemble.Nensembles));

% Show core ensemble neurons
figure;
if(isfield(data,'image'))
    if(size(data.image,3)==3)
        image = data.image(:,:,2);
    else
        image = data.image;
    end
end
imagesc(image); colormap('gray');axis image;
hold on

c = regionprops(data.L,'Centroid');
c = reshape([c.Centroid],2,data.N);
x = c(1,:);
y = c(2,:);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

for i=1:data.N
    plot(x(i),y(i),'ro','MarkerFaceColor','r');
end
for i=1:length(data.NetworkEnsemble.CoreEnsembles)
    ens = data.NetworkEnsemble.CoreEnsembles{i};
    plot(x(ens),y(ens),'go','MarkerFaceColor','g');
end
title(sprintf('Core Ensemble\nNeurons belonging to the core ensemble are colored in green, rest in red'));
% --------------------------------------------------------------------
function controllability_Callback(hObject, eventdata, handles)
% hObject    handle to controllability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function controllability_CC_Callback(hObject, eventdata, handles)
% hObject    handle to controllability_CC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
try
    drivernodes = data.FC.CC.Controllability.drivernodes;
    % Controllability using cross-correlation adjacency matrix
    figure;
    if(isfield(data,'image'))
        if(size(data.image,3)==3)
            image = data.image(:,:,2);
        else
            image = data.image;
        end
    end
    imagesc(image); colormap('gray');axis image;
    hold on
    
    c = regionprops(data.L,'Centroid');
    c = reshape([c.Centroid],2,data.N);
    x = c(1,:);
    y = c(2,:);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    
    for i=1:data.N
        plot(x(i),y(i),'ro','MarkerFaceColor','r');
    end
    
    plot(x(drivernodes),y(drivernodes),'go','MarkerFaceColor','g');
    
    title(sprintf('Network controllability\nDriver nodes colored in green, rest in red\n%d driver nodes, %d unique configurations',data.FC.CC.Controllability.Nd,data.FC.CC.Controllability.Nconfigs));
catch
    errordlg('Error: make sure functional connectivity using cross-correlation exists. Check Preferences and re-compute','Bad input','modal');
    return;
end
% --------------------------------------------------------------------
function controllability_PC_Callback(hObject, eventdata, handles)
% hObject    handle to controllability_PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
try
    drivernodes = data.FC.PC.Controllability.drivernodes;
    % Controllability using partial-correlation adjacency matrix
    figure;
    if(isfield(data,'image'))
        if(size(data.image,3)==3)
            image = data.image(:,:,2);
        else
            image = data.image;
        end
    end
    imagesc(image); colormap('gray');axis image;
    hold on
    
    c = regionprops(data.L,'Centroid');
    c = reshape([c.Centroid],2,data.N);
    x = c(1,:);
    y = c(2,:);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    
    for i=1:data.N
        plot(x(i),y(i),'ro','MarkerFaceColor','r');
    end
    
    plot(x(drivernodes),y(drivernodes),'go','MarkerFaceColor','g');
    
    title(sprintf('Network controllability\nDriver nodes colored in green, rest in red\n%d driver nodes, %d unique configurations',data.FC.PC.Controllability.Nd,data.FC.PC.Controllability.Nconfigs));
catch
    errordlg('Error: make sure functional connectivity using partial-correlation exists. Check Preferences and re-compute','Bad input','modal');
    return;
end
% --------------------------------------------------------------------
function controllability_phase_Callback(hObject, eventdata, handles)
% hObject    handle to controllability_phase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
try
    drivernodes = data.FC.phase.Controllability.drivernodes;
    % Controllability using instantaneous phase adjacency matrix
    figure;
    if(isfield(data,'image'))
        if(size(data.image,3)==3)
            image = data.image(:,:,2);
        else
            image = data.image;
        end
    end
    imagesc(image); colormap('gray');axis image;
    hold on
    
    c = regionprops(data.L,'Centroid');
    c = reshape([c.Centroid],2,data.N);
    x = c(1,:);
    y = c(2,:);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    
    for i=1:data.N
        plot(x(i),y(i),'ro','MarkerFaceColor','r');
    end
    
    plot(x(drivernodes),y(drivernodes),'go','MarkerFaceColor','g');
    
    title(sprintf('Network controllability\nDriver nodes colored in green, rest in red\n%d driver nodes, %d unique configurations',data.FC.phase.Controllability.Nd,data.FC.phase.Controllability.Nconfigs));
catch
    errordlg('Error: make sure functional connectivity using instantaneous phase exists. Check Preferences and re-compute','Bad input','modal');
    return;
end

% --------------------------------------------------------------------
function controllability_GC_Callback(hObject, eventdata, handles)
% hObject    handle to controllability_GC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
try
    drivernodes = data.FC.GC.Controllability.drivernodes;
    % Controllability using Granger causality adjacency matrix
    figure;
    if(isfield(data,'image'))
        if(size(data.image,3)==3)
            image = data.image(:,:,2);
        else
            image = data.image;
        end
    end
    imagesc(image); colormap('gray');axis image;
    hold on
    
    c = regionprops(data.L,'Centroid');
    c = reshape([c.Centroid],2,data.N);
    x = c(1,:);
    y = c(2,:);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    
    for i=1:data.N
        plot(x(i),y(i),'ro','MarkerFaceColor','r');
    end
    
    plot(x(drivernodes),y(drivernodes),'go','MarkerFaceColor','g');
    
    title(sprintf('Network controllability\nDriver nodes colored in green, rest in red\n%d driver nodes, %d unique configurations',data.FC.GC.Controllability.Nd,data.FC.GC.Controllability.Nconfigs));
catch
    errordlg('Error: make sure functional connectivity using Granger causality exists. Check Preferences and re-compute','Bad input','modal');
    return;
end
% --------------------------------------------------------------------
function controllability_TE_Callback(hObject, eventdata, handles)
% hObject    handle to controllability_TE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(handles.curr_folder_idx ==0 || handles.curr_file_idx ==0)
    errordlg('Please first select a folder and a file from the two listboxes.','Bad Input','modal')
    return;
end

try
    load([handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat'],'data');
    handles.data = data;
    
catch
    msg = sprintf('Could not load %s. Please make sure file exists and you used transfer entropy as the method for functional connectivity estimate. Run "Analysis->Recompute"\n',[handles.curr_folder '/analysis-' handles.curr_file(1:end-4) '.mat']);
    errordlg(msg,'Cannot locate file','modal');
    return;
end
try
    drivernodes = data.FC.TE.Controllability.drivernodes;
    % Controllability using transfer entropy adjacency matrix
    figure;
    if(isfield(data,'image'))
        if(size(data.image,3)==3)
            image = data.image(:,:,2);
        else
            image = data.image;
        end
    end
    imagesc(image); colormap('gray');axis image;
    hold on
    
    c = regionprops(data.L,'Centroid');
    c = reshape([c.Centroid],2,data.N);
    x = c(1,:);
    y = c(2,:);
    set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
    
    for i=1:data.N
        plot(x(i),y(i),'ro','MarkerFaceColor','r');
    end
    
    plot(x(drivernodes),y(drivernodes),'go','MarkerFaceColor','g');
    
    title(sprintf('Network controllability\nDriver nodes colored in green, rest in red\n%d driver nodes, %d unique configurations',data.FC.TE.Controllability.Nd,data.FC.TE.Controllability.Nconfigs));
catch
    errordlg('Error: make sure functional connectivity using transfer entropy exists. Check Preferences and re-compute','Bad input','modal');
    return;
end


% --------------------------------------------------------------------
function parallel_computing_Callback(hObject, eventdata, handles)
% hObject    handle to parallel_computing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    choice = questdlg('Do you want to use parallel computing? Requires proper installation of the parallel computing toolbox.','Enable parallel computing?','Enable parallel computing','Disable parallel computing','Enable parallel computing');
    load('params.mat');
    switch choice
        case 'Enable parallel computing'
            params.parallel = 1;
        case 'Disable parallel computing'
            params.parallel = 0;
    end
    save('FluoroSNNAP_code/params.mat','params');
catch
    return;
end


% --------------------------------------------------------------------
function analysis_modules_Callback(hObject, eventdata, handles)
% hObject    handle to analysis_modules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
analysis_options;
