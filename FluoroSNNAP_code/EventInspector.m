function varargout = EventInspector(varargin)
% EVENTINSPECTOR MATLAB code for EventInspector.fig
%      EVENTINSPECTOR, by itself, creates a new EVENTINSPECTOR or raises the existing
%      singleton*.
%
%      H = EVENTINSPECTOR returns the handle to a new EVENTINSPECTOR or the handle to
%      the existing singleton*.
%
%      EVENTINSPECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTINSPECTOR.M with the given input arguments.
%
%      EVENTINSPECTOR('Property','Value',...) creates a new EVENTINSPECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EventInspector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EventInspector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EventInspector

% Last Modified by GUIDE v2.5 23-Dec-2014 12:00:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EventInspector_OpeningFcn, ...
    'gui_OutputFcn',  @EventInspector_OutputFcn, ...
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


% --- Executes just before EventInspector is made visible.
function EventInspector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EventInspector (see VARARGIN)

% Choose default command line output for EventInspector
handles.output = hObject;
set(handles.figure1,'Visible','on');
handles.data = varargin{1};
handles.curr_ROI = [];
imagesc(handles.data.image,'Parent',handles.axes1); colormap('gray');axis(handles.axes1,'image');
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
handles.t = 0:1/handles.data.fps:size(handles.data.F_cell,2)/handles.data.fps-1/handles.data.fps;
axis(handles.axes1,'image');
linkaxes([handles.axes2,handles.axes3,handles.axes6],'x');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EventInspector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EventInspector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function tools_Callback(hObject, eventdata, handles)
% hObject    handle to tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function selectROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=zoom;
setAxesZoomMotion(h, handles.axes2, 'horizontal');
setAxesZoomMotion(h, handles.axes3, 'horizontal');

h = impoint(handles.axes1);
coords = h.getPosition;
handles.curr_ROI = handles.data.L(floor(coords(2)),floor(coords(1)));
delete(h);
if(handles.curr_ROI==0)
    errordlg('ROI not selected. Please try again (try zooming in)','Bad input','modal');
    return;
end
imagesc(handles.data.image,'Parent',handles.axes1); colormap('gray');
hold(handles.axes1,'on');
nROIs = max(handles.data.L(:));
C = regionprops(handles.data.L,'Centroid');
for i=1:nROIs
    B = bwboundaries(handles.data.L==i);
    if(i==handles.curr_ROI)
        for k=1:length(B)
            plot(handles.axes1,B{k}(:,2),B{k}(:,1),'b','LineWidth',3);
        end
    else
        for k=1:length(B)
            plot(handles.axes1,B{k}(:,2),B{k}(:,1),'r');
        end
    end
    text(C(i).Centroid(1),C(i).Centroid(2),num2str(i),'Color','r','Parent',handles.axes1);
end
hold(handles.axes1,'off');
set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);
drawnow;

%Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
%in axes6
plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);

plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
try
    plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
    xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
end
% If overlay spike is checked:
state = get(handles.checkbox1,'Value');
if(state)
    % Overlay spikes
    hold(handles.axes2,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes2);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
    end
    hold(handles.axes2,'off');
    xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
    
    hold(handles.axes6,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes6);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
    end
    hold(handles.axes6,'off');
    xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
     if(handles.data.params.analyze.spike_probability)
        hold(handles.axes3,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes3);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
        end
        hold(handles.axes3,'off');
        xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function enter_number_Callback(hObject, eventdata, handles)
% hObject    handle to enter_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


answer = inputdlg('Enter ROI number to select','Enter ROI',1,{'1'});
handles.curr_ROI = str2num(answer{1});
imagesc(handles.data.image,'Parent',handles.axes1); colormap('gray');
hold(handles.axes1,'on');
nROIs = max(handles.data.L(:));
C = regionprops(handles.data.L,'Centroid');
handles.axes1;
for i=1:nROIs
    B = bwboundaries(handles.data.L==i);
    if(i==handles.curr_ROI)
        for k=1:length(B)
            plot(handles.axes1,B{k}(:,2),B{k}(:,1),'b','LineWidth',3);
        end
    else
        for k=1:length(B)
            plot(handles.axes1,B{k}(:,2),B{k}(:,1),'r');
        end
    end
    text(C(i).Centroid(1),C(i).Centroid(2),num2str(i),'Color','r','Parent',handles.axes1);
end
hold(handles.axes1,'off');
set(handles.axes1,'XTickLabel',[]); set(handles.axes1,'YTickLabel',[]);
drawnow;

%Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
%in axes6
plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);

plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
try
    plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
    xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
end
% If overlay spike is checked:
state = get(handles.checkbox1,'Value');
if(state)
    % Overlay spikes
    hold(handles.axes2,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes2);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
    end
    hold(handles.axes2,'off');
    xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
    
    hold(handles.axes6,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes6);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
    end
    hold(handles.axes6,'off');
    xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
     if(handles.data.params.analyze.spike_probability)
        hold(handles.axes3,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes3);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
        end
        hold(handles.axes3,'off');
        xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end
guidata(hObject,handles);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
state = get(hObject,'Value');
if(isempty(handles.curr_ROI))
    errordlg('You must select an ROI first','Bad input','modal');
    return;
end
if(state)
    % Overlay spikes
    hold(handles.axes2,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes2);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
    end
    hold(handles.axes2,'off');
    xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
    
    hold(handles.axes6,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes6);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
    end
    hold(handles.axes6,'off');
    xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
    if(handles.data.params.analyze.spike_probability)
        hold(handles.axes3,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes3);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
        end
        hold(handles.axes3,'off');
        xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end

if(~state)
    plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);
    
    plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
     if(handles.data.params.analyze.spike_probability)
        plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
        xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end



% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function remove_event_Callback(hObject, eventdata, handles)
% hObject    handle to remove_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[x,~]=getline(handles.axes2);
fr = x.*handles.data.fps;
sp = handles.data.Spikes_cell{handles.curr_ROI};

% Delete any spikes that fall within fr range
idx=find(sp>fr(1) & sp<fr(2));
sp(idx) = [];
handles.data.Spikes_cell{handles.curr_ROI} = sp;
%Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
%in axes6
plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);

plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
try
    plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
    xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
end
% If overlay spike is checked:
state = get(handles.checkbox1,'Value');
if(state)
    % Overlay spikes
    hold(handles.axes2,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes2);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
    end
    hold(handles.axes2,'off');
    xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
    
    hold(handles.axes6,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes6);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
    end
    hold(handles.axes6,'off');
    xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
     if(handles.data.params.analyze.spike_probability)
        hold(handles.axes3,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes3);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
        end
        hold(handles.axes3,'off');
        xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function add_event_Callback(hObject, eventdata, handles)
% hObject    handle to add_event (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,~] = getpts(handles.axes2);
sp = handles.data.Spikes_cell{handles.curr_ROI};
sp = unique([sp ceil(x.*handles.data.fps)']);
handles.data.Spikes_cell{handles.curr_ROI} = sp;
%Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
%in axes6
plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);

plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
try
    plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
    xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
end
% If overlay spike is checked:
state = get(handles.checkbox1,'Value');
if(state)
    % Overlay spikes
    hold(handles.axes2,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes2);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
    end
    hold(handles.axes2,'off');
    xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
    
    hold(handles.axes6,'on');
    spk = handles.data.Spikes_cell{handles.curr_ROI};
    ymn = ylim(handles.axes6);
    for i=1:length(spk)
        line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
    end
    hold(handles.axes6,'off');
    xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
     if(handles.data.params.analyze.spike_probability)
        hold(handles.axes3,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes3);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
        end
        hold(handles.axes3,'off');
        xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
end
guidata(hObject,handles);
% --------------------------------------------------------------------
function new_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to new_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preferences_Callback(hObject, eventdata, handles);
multiWaitbar('CloseAll','Name','');
multiWaitbar('Updating event detection','Busy','Name','Please wait');
spks = EventDetection(handles.data.dF_cell(handles.curr_ROI,:));
handles.data.Spikes_cell{handles.curr_ROI} = spks;
multiWaitbar('CloseAll','Name','');
try
    %Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
    %in axes6
    plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);
    
    plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
    try
        plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
        xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
    % If overlay spike is checked:
    state = get(handles.checkbox1,'Value');
    if(state)
        % Overlay spikes
        hold(handles.axes2,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes2);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
        end
        hold(handles.axes2,'off');
        xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
        
        hold(handles.axes6,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes6);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
        end
        hold(handles.axes6,'off');
        xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
         if(handles.data.params.analyze.spike_probability)
            hold(handles.axes3,'on');
            spk = handles.data.Spikes_cell{handles.curr_ROI};
            ymn = ylim(handles.axes3);
            for i=1:length(spk)
                line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
            end
            hold(handles.axes3,'off');
            xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
        end
    end
end
guidata(hObject,handles);
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

% Save handles.data
prompt = questdlg('Save changes?','Save?','Yes','No','Yes');
switch prompt
    case 'Yes'
        fname = handles.data.filename;
        [folder,file] = fileparts(fname);
        data = handles.data;
        save([folder '/analysis-' file '.mat'],'data');
        msgbox('Analysis updated with new event detections. Please run Analysis -> Recompute in the FluoroSNNAP window to update synchronization and functional connectivity analysis');
        
end
delete(hObject);


% --------------------------------------------------------------------
function new_thr_all_Callback(hObject, eventdata, handles)
% hObject    handle to new_thr_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
preferences_Callback(hObject, eventdata, handles);
multiWaitbar('CloseAll','Name','');
multiWaitbar('Updating event detection for all ROIs.',0,'Name','Please wait');
for i=1:handles.data.N
    spks = EventDetection(handles.data.dF_cell(i,:));
    handles.data.Spikes_cell{i} = spks;
    multiWaitbar('Updating event detection for all ROIs.',i/handles.data.N);
    
end
multiWaitbar('CloseAll','Name','');
try
    %Plot calcium activity in axes2 and inferred spikes in axes3 and DeltaF/F
    %in axes6
    plot(handles.axes2,handles.t,handles.data.F_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes2,'Time (s)'); ylabel(handles.axes2,'Calcium fluorescence (a.u.)'); title(handles.axes2,['Current ROI = ' num2str(handles.curr_ROI)]);
    
    plot(handles.axes6,handles.t,handles.data.dF_cell(handles.curr_ROI,:),'k');
    xlabel(handles.axes6,'Time (s)'); ylabel(handles.axes6,'\DeltaF/F');
    try
        plot(handles.axes3,handles.t,handles.data.foopsi(handles.curr_ROI,:),'k');
        xlabel(handles.axes3,'Time (s)'); ylabel(handles.axes3,'Inferred spikes (a.u.)');
    end
    % If overlay spike is checked:
    state = get(handles.checkbox1,'Value');
    if(state)
        % Overlay spikes
        hold(handles.axes2,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes2);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes2);
        end
        hold(handles.axes2,'off');
        xlabel(handles.axes2,'Time (s)');ylabel(handles.axes2,'Calcium fluorescence (a.u.)');
        
        hold(handles.axes6,'on');
        spk = handles.data.Spikes_cell{handles.curr_ROI};
        ymn = ylim(handles.axes6);
        for i=1:length(spk)
            line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes6);
        end
        hold(handles.axes6,'off');
        xlabel(handles.axes6,'Time (s)');ylabel(handles.axes6,'\DeltaF/F');
         if(handles.data.params.analyze.spike_probability)
            hold(handles.axes3,'on');
            spk = handles.data.Spikes_cell{handles.curr_ROI};
            ymn = ylim(handles.axes3);
            for i=1:length(spk)
                line([handles.t(spk(i)),handles.t(spk(i))],ymn,'Color','b','Parent',handles.axes3);
            end
            hold(handles.axes3,'off');
            xlabel(handles.axes3,'Time (s)');ylabel(handles.axes3,'Inferred spikes (a.u.)');
        end
    end
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function delete_waveform_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to delete_waveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
remove_event_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
add_event_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function preferences_Callback(hObject, eventdata, handles)
% hObject    handle to preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msg1 = sprintf('Event detection method. 1 = template, 2 = deconvolution');
msg2 = sprintf('Threshold for event detection. Enter a number between 0 and 1 for event-based and a scalar >=1 for deconvolution-based method');
msg3 = sprintf('Minimum deltaF/F amplitude to be considered a calcium transient');
load('params.mat');
try
    answer = inputdlg({msg1,msg2,msg3},'Event detection preferences',1,{num2str(params.event_type),num2str(params.event_thresh),num2str(params.event_amplitude)});
    params.event_type=str2num(answer{1});
    params.event_thresh=str2num(answer{2});
    params.event_amplitude=str2num(answer{3});
    save('FluoroSNNAP_code/params.mat','params');
catch
    return;
end
