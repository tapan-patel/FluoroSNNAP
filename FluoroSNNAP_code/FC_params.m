function varargout = FC_params(varargin)
% FC_PARAMS MATLAB code for FC_params.fig
%      FC_PARAMS, by itself, creates a new FC_PARAMS or raises the existing
%      singleton*.
%
%      H = FC_PARAMS returns the handle to a new FC_PARAMS or the handle to
%      the existing singleton*.
%
%      FC_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FC_PARAMS.M with the given input arguments.
%
%      FC_PARAMS('Property','Value',...) creates a new FC_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FC_params_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FC_params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FC_params

% Last Modified by GUIDE v2.5 12-Jan-2015 19:46:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FC_params_OpeningFcn, ...
    'gui_OutputFcn',  @FC_params_OutputFcn, ...
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


% --- Executes just before FC_params is made visible.
function FC_params_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FC_params (see VARARGIN)

% Choose default command line output for FC_params
handles.output = hObject;
load('params.mat');
handles.params = params;

% Load defaults
try
    set(handles.edit1,'String',num2str(handles.params.FC.CC.Nsur));
catch
    handles.params.FC.CC.Nsur = 100;
end
try
    set(handles.edit2,'String',num2str(handles.params.FC.CC.maxlag));
catch
    handles.params.FC.CC.maxlag = .5;
end
try
    set(handles.edit3,'String',num2str(handles.params.FC.PC.alpha));
catch
    
    handles.params.FC.PC.alpha = 0.001;
end
try
    set(handles.edit4,'String',num2str(handles.params.FC.phase.Nsur));
catch
    handles.params.FC.phase.Nsur = 100;
end
try
    set(handles.edit5,'String',num2str(handles.params.FC.phase.alpha));
catch
    handles.params.FC.phase.alpha = 0.001;
end
try
    set(handles.edit6,'String',num2str(handles.params.FC.GC.morder));
catch
    handles.params.FC.GC.morder=20;
end
try
    set(handles.edit7,'String',num2str(handles.params.FC.GC.alpha));
catch
    handles.params.FC.GC.alpha=0.05;
end
handles.params.FC.GC.iter=100;
handles.params.FC.TE.lags=1:10;
try
    set(handles.edit8,'String',num2str(handles.params.FC.TE.Nsur));
catch
    handles.params.FC.TE.Nsur=100;
end
try
    set(handles.popupmenu1,'Value',handles.params.FC.method_idx+1);
catch
    handles.params.FC.method_idx = 6; % Default is to compute FC using all methods
    set(handles.popupmenu1,'Value',handles.params.FC.method_idx+1);
end

try
    switch handles.params.FC.method_idx+1
        case 1
            set(handles.uipanel1,'HighlightColor','w');
            set(handles.uipanel2,'HighlightColor','w');
            set(handles.uipanel3,'HighlightColor','w');
            set(handles.uipanel4,'HighlightColor','w');
            set(handles.uipanel5,'HighlightColor','w');
        case 2
            set(handles.uipanel1,'HighlightColor','r');
        case 3
            set(handles.uipanel2,'HighlightColor','r');
        case 4
            set(handles.uipanel3,'HighlightColor','r');
        case 5
            set(handles.uipanel4,'HighlightColor','r');
        case 6
            set(handles.uipanel5,'HighlightColor','r');
        case 7
            set(handles.uipanel1,'HighlightColor','r');
            set(handles.uipanel2,'HighlightColor','r');
            set(handles.uipanel3,'HighlightColor','r');
            set(handles.uipanel4,'HighlightColor','r');
            set(handles.uipanel5,'HighlightColor','r');
    end
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FC_params wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FC_params_OutputFcn(hObject, eventdata, handles)
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
params = handles.params;
save('FluoroSNNAP_code/params','params');
msgbox('Parameters updated. Exit this GUI to return to FluoroSNNAP');

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

contents = cellstr(get(hObject,'String'));
handles.params.FC.method = contents{get(hObject,'Value')};
handles.params.FC.method_idx = get(hObject,'Value')-1;
switch get(hObject,'Value')
    case 1
        set(handles.uipanel1,'HighlightColor','w');
        set(handles.uipanel2,'HighlightColor','w');
        set(handles.uipanel3,'HighlightColor','w');
        set(handles.uipanel4,'HighlightColor','w');
        set(handles.uipanel5,'HighlightColor','w');
    case 2
        set(handles.uipanel1,'HighlightColor','r');
        set(handles.uipanel2,'HighlightColor','w');
        set(handles.uipanel3,'HighlightColor','w');
        set(handles.uipanel4,'HighlightColor','w');
        set(handles.uipanel5,'HighlightColor','w');
    case 3
        set(handles.uipanel2,'HighlightColor','r');
        set(handles.uipanel3,'HighlightColor','w');
        set(handles.uipanel4,'HighlightColor','w');
        set(handles.uipanel5,'HighlightColor','w');
        set(handles.uipanel1,'HighlightColor','w');
    case 4
        set(handles.uipanel3,'HighlightColor','r');
        set(handles.uipanel4,'HighlightColor','w');
        set(handles.uipanel5,'HighlightColor','w');
        set(handles.uipanel1,'HighlightColor','w');
        set(handles.uipanel2,'HighlightColor','w');
    case 5
        set(handles.uipanel4,'HighlightColor','r');
        set(handles.uipanel5,'HighlightColor','w');
        set(handles.uipanel1,'HighlightColor','w');
        set(handles.uipanel2,'HighlightColor','w');
        set(handles.uipanel3,'HighlightColor','w');
        
    case 6
        set(handles.uipanel5,'HighlightColor','r');
    case 7
        set(handles.uipanel1,'HighlightColor','r');
        set(handles.uipanel2,'HighlightColor','r');
        set(handles.uipanel3,'HighlightColor','r');
        set(handles.uipanel4,'HighlightColor','r');
        set(handles.uipanel5,'HighlightColor','r');
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.params.FC.CC.Nsur = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.params.FC.CC.maxlag = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
handles.params.FC.phase.Nsur = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

handles.params.FC.phase.alpha = str2double(get(hObject,'String'));
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.params.FC.PC.alpha = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
handles.params.FC.GC.morder = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

handles.params.FC.GC.alpha = str2double(get(hObject,'String'));
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
idx = get(hObject,'Value');
switch(idx)
    case 3
        handles.params.FC.GC.iter = -1;
    case 2
        answer= inputdlg('Enter number of iterations per neuron','Iterations',1,{'100'});
        handles.params.FC.GC.iter = str2double(answer{1});
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double
handles.params.FC.TE.Nsur = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double
handles.params.FC.TE.lags = eval(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
