function varargout = ViewImageStack(varargin)
% VIEWIMAGESTACK MATLAB code for ViewImageStack.fig
%      VIEWIMAGESTACK, by itself, creates a new VIEWIMAGESTACK or raises the existing
%      singleton*.
%
%      H = VIEWIMAGESTACK returns the handle to a new VIEWIMAGESTACK or the handle to
%      the existing singleton*.
%
%      VIEWIMAGESTACK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWIMAGESTACK.M with the given input arguments.
%
%      VIEWIMAGESTACK('Property','Value',...) creates a new VIEWIMAGESTACK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewImageStack_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewImageStack_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewImageStack

% Last Modified by GUIDE v2.5 05-Sep-2014 15:33:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ViewImageStack_OpeningFcn, ...
    'gui_OutputFcn',  @ViewImageStack_OutputFcn, ...
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


% --- Executes just before ViewImageStack is made visible.
function ViewImageStack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewImageStack (see VARARGIN)

% Choose default command line output for ViewImageStack
global cancel;
cancel = 0;
handles.output = hObject;
set(handles.figure1,'Visible','on');
if(nargin)
    if(ischar(varargin{1}))
        % This is a file name. Make sure it exists
        filename = varargin{1};
        if(~exist(filename,'file'))
            errordlg([filename ' does not exist. Make sure to include full path to this file.']);
            return
        end
        try
            % Message box
            msg = sprintf('Please wait while reading %s.\nThis may take seconds to minutes',filename);
            h = PleaseWait(msg);
            figure(h);
            [handles.Stack,handles.time,handles.acq_fps] = ReadTiffStack(filename);
            close(h);
        catch
            errordlg(['Could not read ' filename '. Make sure the file is readable through ReadTiffStack.m']);
        end
    else
        % Input is a 3D matrix. Need 3 inputs. If not, set time and acq_fps
        % to Inf.
        handles.Stack = varargin{1};
        try
            handles.time = varargin{2};
        catch
            % assume frame rate is 10
            handles.time = 0:1/10:size(handles.Stack,3)/10-1/10;
        end
        try
            handles.acq_fps = varargin{3};
        catch
            handles.acq_fps = 10;
        end
    end
end
% Show the first image frame in jet colormap. Adjust contrast for all
% images
imin = min(handles.Stack(:));
imax = max(handles.Stack(:));

handles.cmin = imin;
handles.cmax = imax;
handles.bitDepth = class(handles.Stack);
handles.enhance_contrast = 0;
I = handles.Stack(:,:,1);
if(handles.enhance_contrast)
    I = imadjust(I,stretchlim([handles.cmin handles.cmax]));
end
handles.I = I;
imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]); colorbar;
title(['Elapsed time = ' num2str(0) ' (s)']);
set(handles.edit3,'String',num2str(round(handles.cmin)));
set(handles.edit2,'String',num2str(round(handles.cmax)));

handles.N = size(handles.Stack,3);
handles.idx = 1;
handles.fps = 20;
% Update slider
set(handles.slider1,'Min',1,'Max',handles.N,'Value',1);
set(handles.slider1,'SliderStep',[1/handles.N 20/handles.N]);
set(handles.axes1,'XTickLabel',[]);
set(handles.axes1,'YTickLabel',[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewImageStack wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewImageStack_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

idx = round(get(hObject,'Value'));
% Display this image in the appropriate colormap
I = handles.Stack(:,:,idx);
handles.I = I;
if(handles.enhance_contrast)
    I = imadjust(I,stretchlim([handles.cmin handles.cmax]));
end

handles.idx = idx;
imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]);
title(['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
colorbar;
set(handles.axes1,'XTickLabel',[]);
set(handles.axes1,'YTickLabel',[]);
set(handles.text1,'String',['Frame #' num2str(idx)]);
guidata(hObject,handles);
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

fps = str2double(get(hObject,'String'));
if(isnan(fps))
    errordlg(['Display frame rate must be a number > 0']);
    return;
else
    handles.fps = fps;
end
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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Play all images in the stack at frame rate specified, start from current
% index
global cancel;
cancel=0;
fps = handles.fps;

for i=handles.idx:handles.N
    if(~cancel)
        if(handles.enhance_contrast)
            I = imadjust(handles.Stack(:,:,i),stretchlim([handles.cmin handles.cmax]));
        else
            I = handles.Stack(:,:,i);
        end
        imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]);
        title(['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
        colorbar;
        handles.idx = i;
        set(handles.axes1,'XTickLabel',[]);
        set(handles.axes1,'YTickLabel',[]);
        set(handles.slider1,'Value',i);
        set(handles.text1,'String',['Frame #' num2str(i)]);
        pause(1/fps);
    else
        guidata(hObject,handles);
        break;
    end
end
guidata(hObject,handles);

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
handles.cmax = round(str2double(get(hObject,'String')));
if(handles.enhance_contrast)
    I = imadjust(handles.I,stretchlim([handles.cmin handles.cmax]));
else
    I = handles.I;
end
imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]);
title(['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
colorbar;
set(handles.axes1,'XTickLabel',[]);
set(handles.axes1,'YTickLabel',[]);
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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
handles.cmin = round(str2double(get(hObject,'String')));
if(handles.enhance_contrast)
    I = imadjust(handles.I,stretchlim([handles.cmin handles.cmax]));
else
    I = handles.I;
end
imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]);
title(['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
colorbar;
set(handles.axes1,'XTickLabel',[]);
set(handles.axes1,'YTickLabel',[]);
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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
contents = cellstr(get(hObject,'String'));
choice = contents{get(hObject,'Value')};
switch choice
    case 'jet'
        colormap(handles.axes1,'jet');
    case 'gray'
        colormap(handles.axes1,'gray');
    case 'green'
        map = zeros(64,3);
        vals = linspace(double(handles.cmin),double(handles.cmax),64);
        map(:,2) = vals'./double(handles.cmax);
        colormap(handles.axes1,map);
end

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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cancel;
cancel = 1;



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

% If user specifies a different acquisiton frame rate, update the
% handles.time property
dt = 1/str2double(get(hObject,'String'));
if(isnan(dt))
    errordlg('Acquisiton frame rate must be number >0');
    return
else
    handles.time = 0:dt:handles.N*dt-dt;
end
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


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Save the time lapse to movie. User can select start and end frame # and
% file format for saving. The current setting for colormap, enhance
% contrast, color limits will be kept..
message = sprintf('You are about to save the current timelapse data as a movie.\nThe current settings - colormap, display frame rate, contrast enhacement etc. - will be used.\nAre you happy with your choices?');
button = questdlg(message,'Confirm settings','Proceed','Go back','Proceed');
switch button
    case 'Proceed'
        prompt = {'Start frame #:',sprintf('End frame #\n(%d total frames)',handles.N)};
        dlg_title = 'Start/end frame #';
        num_lines = 1;
        def = {'1',num2str(handles.N)};
        try
        answer = inputdlg(prompt,dlg_title,num_lines,def);
       catch
        return;
    end
        start_idx = str2num(answer{1});
        end_idx = str2num(answer{2});
        
        [filename, pathname] = uiputfile(...
            {'*.avi';'*.mp4'},...
            'Save movie as');
        if(end_idx<start_idx || isnan(start_idx) || isnan(end_idx))
            errordlg('Start frame and end frame must be > 0 and end frame cannot be less than start frame.');
        end
        hfig = figure('Visible','off');
        hax = gca;
        map = colormap(handles.axes1);
        colormap(hax,map);
        set(hfig,'PaperPositionMode','auto');
        vidout = VideoWriter(fullfile(pathname, filename));
        vidout.FrameRate = str2double(get(handles.edit1,'String'));
        open(vidout);
        
        hwait = waitbar(0,sprintf('Encoding video\n0%%'),'Name','Please wait. Encoding video',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(hwait,'canceling',0)
        for i=start_idx:end_idx
            if getappdata(hwait,'canceling')
                close(vidout);
            
                break
            end
            I = handles.Stack(:,:,i);
            if(handles.enhance_contrast)
                I = imadjust(I,stretchlim([handles.cmin,handles.cmax]));
            end;
            imagesc(I,'Parent',hax,[handles.cmin handles.cmax]);
            title(hax,['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
            colorbar();
            set(hax,'XTickLabel',[]);
            set(hax,'YTickLabel',[]);
            
            cdata = hardcopy(hfig, '-Dzbuffer', '-r0');
            writeVideo(vidout,cdata);
            waitbar(i/(end_idx-start_idx+1),hwait,sprintf('Encoding video\n%d%%',round(100*i/(end_idx-start_idx+1))));
        end
        close(vidout);
        delete(hfig);
        delete(hwait);

        msgbox(['Video successfully saved to ' fullfile(pathname,filename)]);
    case 'Go back'
        return;
end

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
val = get(hObject,'Value');
if(val)
    % Enhance contrast
    handles.enhance_contrast = 1;
else
    handles.enhance_contrast = 0;
end
if(handles.enhance_contrast)
    I = imadjust(handles.I,stretchlim([handles.cmin handles.cmax]));
else
    I = handles.I;
end
imagesc(I,'Parent',handles.axes1,[handles.cmin handles.cmax]);
title(['Elapsed time = ' num2str(handles.time(handles.idx)) ' (s)']);
colorbar;
set(handles.axes1,'XTickLabel',[]);
set(handles.axes1,'YTickLabel',[]);
guidata(hObject,handles);
