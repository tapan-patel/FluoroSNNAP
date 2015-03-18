function varargout = analysis_options(varargin)
% ANALYSIS_OPTIONS MATLAB code for analysis_options.fig
%      ANALYSIS_OPTIONS, by itself, creates a new ANALYSIS_OPTIONS or raises the existing
%      singleton*.
%
%      H = ANALYSIS_OPTIONS returns the handle to a new ANALYSIS_OPTIONS or the handle to
%      the existing singleton*.
%
%      ANALYSIS_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSIS_OPTIONS.M with the given input arguments.
%
%      ANALYSIS_OPTIONS('Property','Value',...) creates a new ANALYSIS_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analysis_options_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analysis_options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analysis_options

% Last Modified by GUIDE v2.5 15-Mar-2015 17:56:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @analysis_options_OpeningFcn, ...
    'gui_OutputFcn',  @analysis_options_OutputFcn, ...
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


% --- Executes just before analysis_options is made visible.
function analysis_options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analysis_options (see VARARGIN)

% Choose default command line output for analysis_options
handles.output = hObject;
set(handles.figure1,'Visible','on');
load('params.mat');
try
    if(params.analyze.deltaF)
        set(handles.checkbox1,'Value',1);
    else
        set(handles.checkbox1,'Value',0);
    end
catch
    params.analyze.deltaF=1;
end
try
    if(params.analyze.detect_events)
        set(handles.checkbox2,'Value',1);
    else
        set(handles.checkbox2,'Value',0);
    end
catch
    params.analyze.detect_events=1;
end
try
    if(params.analyze.sca)
        set(handles.checkbox3,'Value',1);
    else
        set(handles.checkbox3,'Value',0);
    end
catch
    params.analyze.sca=1;
end
try
    if(params.analyze.FC)
        set(handles.checkbox4,'Value',1);
    else
        set(handles.checkbox4,'Value',0);
    end
catch
    params.analyze.FC=1;
end

try
    if(params.analyze.controllability)
        set(handles.checkbox5,'Value',1);
    else
        set(handles.checkbox5,'Value',0);
    end
catch
    params.analyze.controllability=1;
end

try
    if(params.analyze.kinetics)
        set(handles.checkbox6,'Value',1);
    else
        set(handles.checkbox6,'Value',0);
    end
catch
    params.analyze.kinetics=1;
end
try
    if(params.analyze.spike_probability)
        set(handles.checkbox7,'Value',1);
    else
        set(handles.checkbox7,'Value',0);
    end
catch
    params.analyze.spike_probability=1;
end
try
    if(params.analyze.ensembles)
        set(handles.checkbox8,'Value',1);
    else
        set(handles.checkbox8,'Value',0);
    end
catch
    params.analyze.ensembles=1;
end
try
    if(params.analyze.figure)
        set(handles.checkbox9,'Value',1);
    else
        set(handles.checkbox9,'Value',0);
    end
catch
    params.analyze.figure=1;
end

handles.params = params;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analysis_options wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = analysis_options_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if(get(hObject,'Value'))
    handles.params.analyze.deltaF = 1;
else
    handles.params.analyze.deltaF = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if(get(hObject,'Value'))
    handles.params.analyze.detect_events = 1;
else
    handles.params.analyze.detect_events = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
if(get(hObject,'Value'))
    handles.params.analyze.sca = 1;
else
    handles.params.analyze.sca = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
if(get(hObject,'Value'))
    handles.params.analyze.FC = 1;
else
    handles.params.analyze.FC = 0;
    handles.params.analyze.controllability = 0;
    
    handles.params.analyze.figure=0;
    set(handles.checkbox5,'Value',0);
  
    set(handles.checkbox9,'Value',0);
end
guidata(hObject,handles);

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5
if(get(hObject,'Value'))
    if(handles.params.analyze.FC==0)
        handles.params.analyze.controllability = 0;
        set(hObject,'Value',0);
        errordlg('You must choose to infer functional connectivity in order to use this module','Bad input','modal');
        return
    else
        handles.params.analyze.controllability = 1;
    end
    
else
    handles.params.analyze.controllability = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6
if(get(hObject,'Value'))
    handles.params.analyze.kinetics = 1;
else
    handles.params.analyze.kinetics = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7
if(get(hObject,'Value'))
    handles.params.analyze.spike_probability = 1;
else
    handles.params.analyze.spike_probability = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8
if(get(hObject,'Value'))
  
        handles.params.analyze.ensembles = 1;
   
    
else
    handles.params.analyze.ensembles = 0;
end
guidata(hObject,handles);

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9

if(get(hObject,'Value'))
    if(handles.params.analyze.FC==0)
        handles.params.analyze.figure = 0;
        set(hObject,'Value',0);
        errordlg('You must choose to infer functional connectivity in order to use this module','Bad input','modal');
        return
    else
        handles.params.analyze.figure = 1;
    end
else
    handles.params.analyze.figure = 0;
end
guidata(hObject,handles);
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save params
params = handles.params;
save('FluoroSNNAP_code/params.mat','params');
delete(handles.figure1);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure1);
