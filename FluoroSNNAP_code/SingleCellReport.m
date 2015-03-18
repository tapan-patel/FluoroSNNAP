function varargout = SingleCellReport(varargin)
% SINGLECELLREPORT MATLAB code for SingleCellReport.fig
%      SINGLECELLREPORT, by itself, creates a new SINGLECELLREPORT or raises the existing
%      singleton*.
%
%      H = SINGLECELLREPORT returns the handle to a new SINGLECELLREPORT or the handle to
%      the existing singleton*.
%
%      SINGLECELLREPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLECELLREPORT.M with the given input arguments.
%
%      SINGLECELLREPORT('Property','Value',...) creates a new SINGLECELLREPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SingleCellReport_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SingleCellReport_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SingleCellReport

% Last Modified by GUIDE v2.5 26-Nov-2014 08:34:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SingleCellReport_OpeningFcn, ...
                   'gui_OutputFcn',  @SingleCellReport_OutputFcn, ...
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


% --- Executes just before SingleCellReport is made visible.
function SingleCellReport_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SingleCellReport (see VARARGIN)

% Choose default command line output for SingleCellReport
handles.output = hObject;
set(handles.figure1,'Visible','on');
handles.data = varargin{1};
handles.curr_ROI = [];
imagesc(handles.data.image,'Parent',handles.axes1); colormap('jet');
hold(handles.axes1,'on');
nROIs = max(handles.data.L(:));
C = regionprops(handles.data.L,'Centroid');
for i=1:nROIs
    B = bwboundaries(handles.data.L==i);
    for k=1:length(B)
        plot(handles.axes1,B{k}(:,2),B{k}(:,1),'r');
    end
   text(C(i).Centroid(1),C(i).Centroid(2),num2str(i),'Color','r','Parent',handles.axes1); 
end
hold(handles.axes1,'off');
set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);


dat = zeros(handles.data.N,10);

params = handles.data.params;
if(params.FC.method_idx==1)
    A = handles.data.FC.CC.A;
    Ci = handles.data.FC.CC.modularity_Ci;
    Q = handles.data.FC.CC.modularity_Q;
    clu = handles.data.FC.CC.clustering_coef;
elseif(params.FC.method_idx==2)
    A = handles.data.FC.PC.A;
    Ci = handles.data.FC.PC.modularity_Ci;
    Q = handles.data.FC.PC.modularity_Q;
    clu = handles.data.FC.PC.clustering_coef;
elseif(params.FC.method_idx==3)
    A = handles.data.FC.phase.A;
    Ci = handles.data.FC.phase.modularity_Ci;
    Q = handles.data.FC.phase.modularity_Q;
    clu = handles.data.FC.phase.clustering_coef;
elseif(params.FC.method_idx==4)
    A = handles.data.FC.GC.A;
    Ci = handles.data.FC.GC.modularity_Ci;
    Q = handles.data.FC.GC.modularity_Q;
    clu = handles.data.FC.GC.clustering_coef;
elseif(params.FC.method_idx==5)
    A = handles.data.FC.TE.peakTE;
    Ci = handles.data.FC.TE.modularity_Ci;
    Q = handles.data.FC.TE.modularity_Q;
    clu = handles.data.FC.TE.clustering_coef;
else
    A = handles.data.FC.phase.A;
    Ci = handles.data.FC.phase.modularity_Ci;
    Q = handles.data.FC.phase.modularity_Q;
    clu = handles.data.FC.phase.clustering_coef;
end
for i=1:handles.data.N
    x = sort(handles.data.F_cell(i,:));
    dat(i,1) = mean(x(1:ceil(.1*length(x))));
    dat(i,2) = numel(handles.data.Spikes_cell{i});
    dat(i,3) = mean(diff(handles.data.Spikes_cell{i}))/handles.data.fps;
    dat(i,4) = handles.data.DF(i);
    dat(i,5) = handles.data.CV(i);
    r = handles.data.rise_time{i}; tau = handles.data.fall_time{i};
    dat(i,6) = mean(r(~isnan(r)));
    dat(i,7) = mean(tau(~isnan(tau)));
    dat(i,8) = sum(A(i,:))/handles.data.N;
    dat(i,9) = max(handles.data.SynchroCluster.PI(i,:));
    dat(i,10) = nnz(handles.data.SynchroCluster.PI(i,:)>.01);
    dat(i,11) = clu(i);
end
handles.dat = dat;
columnname = {'Baseline fluorescence','Total events','<ISI> (s)','Amplitude','CV','Rise time (s)','Fall time (s)','Connectivity index','Participation index','# assemblies','Clustering coef'};
columneditable = [false false false false false false false false false false false];

set(handles.uitable1,'Data', dat,'ColumnName', columnname,'ColumnEditable', columneditable);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SingleCellReport wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SingleCellReport_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
