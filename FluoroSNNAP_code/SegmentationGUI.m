function varargout = SegmentationGUI(varargin)
% SEGMENTATIONGUI MATLAB code for SegmentationGUI.fig
%      SEGMENTATIONGUI, by itself, creates a new SEGMENTATIONGUI or raises the existing
%      singleton*.
%
%      H = SEGMENTATIONGUI returns the handle to a new SEGMENTATIONGUI or the handle to
%      the existing singleton*.
%
%      SEGMENTATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTATIONGUI.M with the given input arguments.
%
%      SEGMENTATIONGUI('Property','Value',...) creates a new SEGMENTATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SegmentationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SegmentationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SegmentationGUI

% Last Modified by GUIDE v2.5 22-Dec-2014 19:17:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SegmentationGUI_OpeningFcn, ...
    'gui_OutputFcn',  @SegmentationGUI_OutputFcn, ...
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


% --- Executes just before SegmentationGUI is made visible.
function SegmentationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SegmentationGUI (see VARARGIN)

set(handles.figure1,'Visible','on');
handles.I = uint16(varargin{1});
handles.filename = varargin{2};
min_max = stretchlim(handles.I)*2^16;
imagesc(handles.I,'Parent',handles.axes1,[min_max(1),min_max(2)]); colormap('gray');
freezeColors;
set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);
handles.L = false(size(handles.I));
handles.Lold = handles.L;
imagesc(handles.L,'Parent',handles.axes2); colormap('gray');
set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
linkaxes([handles.axes1,handles.axes2],'xy');
handles.curr_ROI = 0;
handles.radius = 20;
% Choose default command line output for SegmentationGUI
handles.output = handles.L;
% set(handles.figure1,'toolbar','figure');
% set(handles.figure1,'menubar','figure');
global stop;
stop=0;
handles.ica = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SegmentationGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SegmentationGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop;
imshow((handles.I),'Parent',handles.axes1);
hold on
while(1)
    h = impoint(handles.axes1);
    mask=false(size(handles.I));
    coord = floor(h.getPosition);
    mask(coord(2)-4:coord(2)+4,coord(1)-4:coord(1)+4)=1;
    bw = activecontour(handles.I,mask);
    handles.Lold = handles.L;
    if(~handles.ica)
        
        handles.L = logical(handles.L) + bw;
        guidata(hObject,handles);
    else
        handles.L(bw) = max(handles.L(:))+1;
        guidata(hObject,handles);
        
    end
    delete(h);
    B = bwboundaries(bw);
    for i=1:length(B)
        plot(B{i}(:,2),B{i}(:,1),'r');
    end
end
guidata(hObject,handles);
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,folder] = uigetfile('','Select a segmentation file to load as either .mat or .tif format');
[~,~,ext] = fileparts(file);
switch ext
    case '.mat'
        data = load(fullfile(folder,file));
        handles.L = data.L;
        if(isfield(data,'ica') && ~isempty(data.ica))
            handles.ica = data.ica;
        end
    case '.tif'
        handles.L = logical(imread(fullfile(folder,file)));
end

guidata(hObject,handles);
Update(handles);

% --------------------------------------------------------------------
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[folder,file] = fileparts(handles.filename);
savename = fullfile(folder,['Segmentation-' file]);
if(~handles.ica)
    L = bwlabel(handles.L);
else
    % Make sure to renumber ROIs so they are contiguous. If for example,
    % the user deletes ROI #50, the rest of ROIs need be renumberd - 51
    % becomes 50, 52 becomes 51 etc.
    for k=1:3
        nROIs = max(handles.L(:));
        for j=1:nROIs
            if(isempty(find(handles.L==j,1))) % This ROI does not exist
                for i=j:nROIs-1
                    handles.L(handles.L==i+1) = i;
                end
                nROIs = nROIs-1;
            end
        end
    end
    L = handles.L;
end
ica = handles.ica;
save(savename,'L','ica');
imwrite(L,[savename '.tif']);
msgbox(['Segmentation file saved to ' savename],'Save successful');
disp(['Segmentation file saved to ' savename]);
guidata(hObject,handles);
% --------------------------------------------------------------------
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Delete segmentations
handles.L = false(size(handles.I));
handles.ica = 0;
Update(handles);
guidata(hObject,handles);

function Update(handles)
min_max = stretchlim(handles.I)*2^16;
imagesc(handles.I,'Parent',handles.axes1,[min_max(1),min_max(2)]);

set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);
hold(handles.axes1,'on');
if(~handles.ica)
    handles.L = bwlabel(handles.L);
    B = bwboundaries(handles.L);
    for i=1:length(B)
        plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
    end
else
    % Make sure to renumber ROIs so they are contiguous. If for example,
    % the user deletes ROI #50, the rest of ROIs need be renumberd - 51
    % becomes 50, 52 becomes 51 etc.
    nROIs = max(handles.L(:));
    for j=1:nROIs
        if(isempty(find(handles.L==j,1))) % This ROI does not exist
            for i=j:nROIs-1
                handles.L(handles.L==i+1) = i;
            end
            nROIs = nROIs-1;
        end
    end
    nROIs = max(handles.L(:));
    for j=1:nROIs
        B = bwboundaries(handles.L==j);
        for i=1:length(B)
            plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
        end
    end
end
hold(handles.axes1,'off')
imshow(label2rgb(handles.L),'Parent',handles.axes2);
title(handles.axes2,sprintf('%d ROIs identified',max(handles.L(:))));

hold(handles.axes2,'on');
C = regionprops(handles.L,'Centroid');
for i=1:length(C)
    text(C(i).Centroid(1),C(i).Centroid(2),num2str(i),'Color','k','FontSize',14,'Parent',handles.axes2);
end
hold(handles.axes2,'off');
% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = questdlg('Exiting. Did you save your work?','Quit?','yes','no','yes');
switch answer
    case 'yes'
        delete(handles.figure1);
    case 'no'
        uicontrol(handles.figure1);
end

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get parameters to threshold
prompt = {'Threshold value:','Minimum ROI area (px)','Maximum ROI area (px)'};
dlg_title = 'Input parameters';
num_lines = 1;
thresh = graythresh(handles.I)*2^16;
def = {num2str(thresh),'50','500'};
try
answer = inputdlg(prompt,dlg_title,num_lines,def);
catch
    return;
end
thresh = str2double(answer{1});
min_area = str2double(answer{2});
max_area = str2double(answer{3});

BW = handles.I>thresh;
imerode(BW,strel('disk',3));
imfill(BW,'holes');
L = bwlabel(BW);
C = regionprops(L,'Area');
A = [C.Area];
idx = find(A>max_area | A<min_area);
for i=1:length(idx)
    BW(L==idx(i))=0;
end
L = handles.L + bwlabel(BW);
handles.Lold = handles.L;
if(~handles.ica)
    handles.L = bwlabel(L);
else
    L = bwlabel(L);
    nROIs = max(L(:));
    for i=1:nROIs
        handles.L(L==i) = max(handles.L(:))+1;
    end
end
guidata(hObject,handles);
Update(handles);
% --------------------------------------------------------------------
function activecontour_Callback(hObject, eventdata, handles)
% hObject    handle to activecontour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ica_Callback(hObject, eventdata, handles)
% hObject    handle to ica (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Automated segmentation using ICA
% This is a GUI wrapper around CellSort code developed by Mukamel,
% publication Neuron 63:747 (2009).

% if(exist('cellsort_preprocessed_data','dir'))
%     !rm -rf cellsort_preprocessed_data
% end
try
prompt = {'Frames to process:','Number of principle components:', 'Number of PCs to use:','Weight of temporal information (0=pure spatial ICA, 1=pure temporal ICA)',...
    'Number of independent components','Standard deviation of Gaussian smoothing kernel (px)','Threshold for spatial filters (s.d.)'...
    'Min and max area of ROIs (px)','Downsample, scalar [temporal spatial]:','Visualize processing (slows computation)'};
dlg_title = 'Input parameters for ICA-based segmentation';
num_lines = 1;
info = imfinfo(handles.filename);
def = {sprintf('[1 %d]',length(info)),'150','1:100','1','100','4','2','[50 500]','[1 1]','0'};

answer = inputdlg(prompt,dlg_title,num_lines,def);
fn = handles.filename;
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

multiWaitbar('CloseAll','Name','');
multiWaitbar('Performing PCA',0,'Name',fn,'CancelFcn', @(a,b) disp( ['Cancel ',a] ) );
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
Centroid = cell(nsegments,1);
cntr = 1;
for i=1:nsegments
    x = logical(squeeze(ica_segments(i,:,:)));
    C = regionprops(x,'Centroid');
    if(~isempty(C))
        
    Lica_stack{cntr} = x;
    Centroid{cntr} = C.Centroid;
    cntr = cntr+1;
    end
end
Lica_stack(cntr:nsegments) = [];
Centroid(cntr:nsegments) = [];
nsegments=length(Lica_stack);
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
        C1 = Centroid{i};
        C2 = Centroid{j};
        D = sqrt( (C1(1)-C2(1))^2 + (C1(2)-C2(2))^2);
        if(D<10)
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
end
remove = unique(remove);
Lica_stack(remove) = [];

nsegments = size(Lica_stack,1);
Lica = zeros(size(handles.I));

for i=1:nsegments
    I = logical(imresize(Lica_stack{i},dsamp(2)));
    Lica(I) = i;
end
handles.Lold = handles.L;
handles.L = Lica;
handles.Lica_stack = Lica_stack;
handles.ica = 1;
guidata(hObject,handles);
Update(handles);
multiWaitbar('Performing segmentation',1,'Color','r');
multiWaitbar('CloseAll');


% --------------------------------------------------------------------
function select_seed_Callback(hObject, eventdata, handles)
% hObject    handle to select_seed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global stop;
stop = 0;
% min_max = stretchlim(handles.I)*2^16;
% imagesc(handles.I,'Parent',handles.axes1,[min_max(1),min_max(2)]); colormap('gray');
% set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);
hold(handles.axes1,'on');
while(1)
    if(stop)
        if(~handles.ica)
            handles.L = bwlabel(handles.L);
        end
        guidata(hObject,handles);
        return;
    else
        
        h = impoint(handles.axes1);
        if(~stop)
            
            mask=false(size(handles.I));
            coord = floor(h.getPosition);
            mask(coord(2)-4:coord(2)+4,coord(1)-4:coord(1)+4)=1;
            bw = activecontour(handles.I,mask);
            if(~handles.ica)
                handles.Lold = bwlabel(handles.L);
                handles.L = logical(handles.L) + bw;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc(bwlabel(handles.L),'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
            else
                handles.Lold = handles.L;
                handles.L(bw) = max(handles.L(:))+1;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc(handles.L,'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
            end
            guidata(hObject,handles);
        else
            delete(h);
        end
    end
end
if(~handles.ica)
    handles.L = bwlabel(handles.L);
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop;
stop = 1;
cla(handles.axes1);
Update(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function edit_Callback(hObject, eventdata, handles)
% hObject    handle to edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function select_Callback(hObject, eventdata, handles)
% hObject    handle to select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(~handles.ica)
    handles.L = bwlabel(handles.L);
    
end
h = impoint(handles.axes1);
coords = floor(h.getPosition);
curr_ROI = handles.L(coords(2),coords(1));

delete(h);
if(curr_ROI~=0)
    % Highlight this ROI so it can be editable. If the ROI is an ellipse,
    % create an ellipse object; otherwise create a polygon object
    C = regionprops(handles.L==curr_ROI,'MajorAxisLength','MinorAxisLength');
    if(C.MajorAxisLength/C.MinorAxisLength>0.9 && C.MajorAxisLength/C.MinorAxisLength<1.1)
        radius = C.MajorAxisLength;
        h = imellipse(handles.axes1,[coords(1)-radius/2 coords(2)-radius/2 radius radius]);
        
    else
        
        B = bwboundaries(handles.L==curr_ROI);
        h = impoly(handles.axes1,[B{1}(1:10:end,2),B{1}(1:10:end,1)]);
        
    end
    handles.h = h;
    handles.curr_ROI = curr_ROI;
    [yrange, xrange] = size(handles.I);
    set(handles.axes1,'XLim',[max([.5 coords(1)-.5*coords(1)]) min([xrange+.5 coords(1)+.5*coords(1)])]);
    zoomylim1 = max([.5 coords(2)-.5*coords(2)]);
    zoomylim2 = min([yrange zoomylim1+range(get(handles.axes1,'XLim'))*yrange/xrange]);
    
    set(handles.axes1,'YLim',[zoomylim1 zoomylim2]);
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function delete_Callback(hObject, eventdata, handles)
% hObject    handle to delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Lold = handles.L;
handles.L(handles.L==handles.curr_ROI) = 0;
guidata(hObject,handles);
Update(handles);

% --------------------------------------------------------------------
function manual_Callback(hObject, eventdata, handles)
% hObject    handle to manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function set_radius_Callback(hObject, eventdata, handles)
% hObject    handle to set_radius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Prompt the user to set the radius of a circle to automatically draw ROIs
try
answer = inputdlg('Enter radius of circle (px):','Radius',1,{'20'});
catch
    return;
end
handles.radius = str2num(answer{1});
guidata(hObject,handles);
% --------------------------------------------------------------------
function freehand_Callback(hObject, eventdata, handles)
% hObject    handle to freehand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global stop;
stop = 0;
hold(handles.axes1,'on');
msg = sprintf('Freehand draw an ROI.\nKeep drawing ROIs.\nClick Tools->Manual->Stop to end selection');
hmsg = msgbox(msg);
uiwait(hmsg);
while(1)
    if(stop)
        break;
    else
       try
        if(~stop)
            
           
            h = imfreehand(handles.axes1);
          
            bw = h.createMask;
            C = regionprops(bw,'Centroid');
            coord = C.Centroid;
            [yrange, xrange] = size(handles.I);
                set(handles.axes1,'XLim',[max([.5 coord(1)-.5*coord(1)]) min([xrange+.5 coord(1)+.5*coord(1)])]);
                zoomylim1 = max([.5 coord(2)-.5*coord(2)]);
                zoomylim2 = min([yrange zoomylim1+range(get(handles.axes1,'XLim'))*yrange/xrange]);
                set(handles.axes1,'YLim',[zoomylim1 zoomylim2]);
            if(~handles.ica)
                handles.Lold = handles.L;
                handles.L = logical(handles.L) + bw;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc(bwlabel(handles.L),'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
                guidata(hObject,handles);
            else
                handles.Lold = handles.L;
                handles.L(bw) = max(handles.L(:))+1;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc(handles.L,'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
                guidata(hObject,handles);
            end
          
        else
            delete(h);
        end
       end
    end
end
if(~handles.ica)
    
    handles.L = bwlabel(handles.L);
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function ellipse_Callback(hObject, eventdata, handles)
% hObject    handle to ellipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global stop;
stop = 0;
radius = handles.radius;
hold(handles.axes1,'on');
msg = sprintf('Click on the center of a cell body. This will place an ellipse in blue.\nResize the ellipse to make sure the ROI is completely deliniated\nDouble-click to accept and move onto the next ROI.\nClick Tools->Manual->Stop to end selection');
h=msgbox(msg);
uiwait(h);
while(1)
    if(stop)
        break;
    else
        h = impoint(handles.axes1);
        if(~stop)
            
            coord = floor(h.getPosition);
            delete(h);
              [yrange, xrange] = size(handles.I);
                set(handles.axes1,'XLim',[max([.5 coord(1)-.5*coord(1)]) min([xrange+.5 coord(1)+.5*coord(1)])]);
                zoomylim1 = max([.5 coord(2)-.5*coord(2)]);
                zoomylim2 = min([yrange zoomylim1+range(get(handles.axes1,'XLim'))*yrange/xrange]);
                set(handles.axes1,'YLim',[zoomylim1 zoomylim2]);
            h = imellipse(handles.axes1,[coord(1)-radius/2 coord(2)-radius/2 radius radius]);
            wait(h);
            bw = h.createMask;
            if(~handles.ica)
                handles.Lold = handles.L;
                handles.L = logical(handles.L) + bw;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc(bwlabel(handles.L),'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
                guidata(hObject,handles);
            else
                handles.Lold = handles.L;
                handles.L(bw) = max(handles.L(:))+1;
                delete(h);
                B = bwboundaries(bw);
                for i=1:length(B)
                    plot(handles.axes1,B{i}(:,2),B{i}(:,1),'r');
                end
                imagesc((handles.L),'Parent',handles.axes2); colormap('gray');
                set(handles.axes2,'XTickLabel',[]); set(handles.axes2,'YTickLabel',[]);
                guidata(hObject,handles);
            end
        else
            delete(h);
        end
    end
end
if(~handles.ica)
    handles.L = bwlabel(handles.L);
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
select_Callback(hObject,eventdata,guidata(hObject))

% --------------------------------------------------------------------
function uipushtool4_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete_Callback(hObject, eventdata, handles)

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('Exiting. Did you save your work?','Quit?','yes','no','yes');
switch answer
    case 'yes'
        delete(handles.figure1);
    case 'no'
        uicontrol(handles.figure1);
end


% --------------------------------------------------------------------
function updateROI_Callback(hObject, eventdata, handles)
% hObject    handle to updateROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Delete the mask of curr_ROI from the label matrix and update it with the
% mask of h
try
handles.L(handles.L==handles.curr_ROI) = 0;
mask = handles.h.createMask;
handles.L(mask) = handles.curr_ROI;
guidata(hObject,handles);
Update(handles);
end

% --------------------------------------------------------------------
function uipushtool5_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Lold = handles.L;
h = imfreehand(handles.axes1);
mask = h.createMask;
handles.L(mask) = 0;
delete(h);

guidata(hObject,handles);
Update(handles);
% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Freehand delete. Useful if there are erroneous processes that are
% segmented together with the cell body and you want to eliminate the
% processes without deleting the entire ROI
handles.Lold = handles.L;
h = imfreehand(handles.axes1);
mask = h.createMask;
handles.L(mask) = 0;
delete(h);
guidata(hObject,handles);
Update(handles);


% --------------------------------------------------------------------
function uipushtool6_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear_Callback(hObject,eventdata,handles);


% --------------------------------------------------------------------
function uipushtool7_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load_Callback(hObject,eventdata,handles);


% --------------------------------------------------------------------
function uipushtool8_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_Callback(hObject,eventdata,handles);


% --------------------------------------------------------------------
function Undo_Callback(hObject, eventdata, handles)
% hObject    handle to Undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.L = handles.Lold;
Update(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function remove_neurites_Callback(hObject, eventdata, handles)
% hObject    handle to remove_neurites (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Cell bodies are round whereas neurites are more straight. For each
% labeled ROI, determine a circularity factor and prompt the user to delete
% ROIs that fall below a certain roundness. Here, 1 = perfect circle and 0
% = line. Circularity = 4*pi*area/perimeter^2

% Prompt user
handles.Lold = handles.L;
try
answer = inputdlg('Circularity threshold (0 = line, 1 = circle)',...
    'Choose circularity threshold',1,{'0.2'});
catch
    return;
end

thresh = str2num(answer{1});
if(~handles.ica)
    C = regionprops(handles.L,'Area','Perimeter');
    Area = [C.Area];
    Perimeter = [C.Perimeter];
    Circularity = 4*pi*Area./(Perimeter.^2);
    idx = find(Circularity<thresh);
    for i=1:length(idx)
        handles.L(handles.L==idx(i)) = 0;
    end
    handles.L = bwlabel(handles.L);
    Update(handles);
else
    % Make sure to renumber ROIs so they are contiguous. If for example,
    % the user deletes ROI #50, the rest of ROIs need be renumberd - 51
    % becomes 50, 52 becomes 51 etc.
    for k=1:3
        nROIs = max(handles.L(:));
        for j=1:nROIs
            if(isempty(find(handles.L==j,1))) % This ROI does not exist
                for i=j:nROIs-1
                    handles.L(handles.L==i+1) = i;
                end
                nROIs = nROIs-1;
            end
        end
    end
    nROIs = max(handles.L(:));
    for i=1:nROIs
        C = regionprops(handles.L==i,'Area','Perimeter');
        A = C.Area; P = C.Perimeter;
        f = 4*pi*A/(P^2);
        if(f<thresh)
            handles.L(handles.L==i) = 0;
        end
    end
    
    nROIs = max(handles.L(:));
    for j=1:nROIs
        if(isempty(find(handles.L==j,1))) % This ROI does not exist
            for i=j:nROIs-1
                handles.L(handles.L==i+1) = i;
            end
            nROIs = nROIs-1;
        end
    end
end
Update(handles);
guidata(hObject,handles);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Resets display magnification
[yrange xrange] = size(handles.I);
set(handles.axes1,'XLim',[.5 xrange+.5]);
set(handles.axes1,'YLim',[.5 yrange+.5]);


% --------------------------------------------------------------------
function shrink_Callback(hObject, eventdata, handles)
% hObject    handle to shrink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    answer = inputdlg('Enter radius (px) to shrink the ROIs by: ','Image erosion radius',1,{'3'});
    radius = str2num(answer{1});
catch
    return;
end

handles.Lold = handles.L;
handles.L = imerode(handles.L,strel('disk',radius));
guidata(hObject,handles);
Update(handles);
