function varargout = gui_getChanLocs(varargin)
% GUI_GETCHANLOCS MATLAB code for gui_getChanLocs.fig
%      GUI_GETCHANLOCS, by itself, creates a new GUI_GETCHANLOCS or raises the existing
%      singleton*.
%
%      H = GUI_GETCHANLOCS returns the handle to a new GUI_GETCHANLOCS or the handle to
%      the existing singleton*.
%
%      GUI_GETCHANLOCS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GETCHANLOCS.M with the given input arguments.
%
%      GUI_GETCHANLOCS('Property','Value',...) creates a new GUI_GETCHANLOCS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_getChanLocs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_getChanLocs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_getChanLocs

% Last Modified by GUIDE v2.5 23-Jan-2018 17:24:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_getChanLocs_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_getChanLocs_OutputFcn, ...
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


% --- Executes just before gui_getChanLocs is made visible.
function gui_getChanLocs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_getChanLocs (see VARARGIN)

% Choose default command line output for gui_getChanLocs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_getChanLocs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_getChanLocs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in anonymizeFace.
function anonymizeFace_Callback(hObject, eventdata, handles)
% hObject    handle to anonymizeFace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of anonymizeFace


% --- Executes on button press in deleteTxtOutput.
function deleteTxtOutput_Callback(hObject, eventdata, handles)
% hObject    handle to deleteTxtOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deleteTxtOutput



function objPath_Callback(hObject, eventdata, handles)
% hObject    handle to objPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of objPath as text
%        str2double(get(hObject,'String')) returns contents of objPath as a double


% --- Executes during object creation, after setting all properties.
function objPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to objPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveName_Callback(hObject, eventdata, handles)
% hObject    handle to saveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveName as text
%        str2double(get(hObject,'String')) returns contents of saveName as a double


% --- Executes during object creation, after setting all properties.
function saveName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browseInput.
function browseInput_Callback(hObject, eventdata, handles)
% hObject    handle to browseInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetdir(pwd);
set(handles.objPath,'String', file)
set(handles.saveName,'String',strcat(file,filesep,'getChanLocs.txt'))


% --- Executes on button press in browseOutput.
function browseOutput_Callback(hObject, eventdata, handles)
% hObject    handle to browseOutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uiputfile('*.txt');
set(handles.saveName,'String', strcat(path,file))



function moveElecInwards_Callback(hObject, eventdata, handles)
% hObject    handle to moveElecInwards (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of moveElecInwards as text
%        str2double(get(hObject,'String')) returns contents of moveElecInwards as a double


% --- Executes during object creation, after setting all properties.
function moveElecInwards_CreateFcn(hObject, eventdata, handles)
% hObject    handle to moveElecInwards (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in help.
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
help getChanLocs

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
% hObject    handle to ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
objPath  = get(handles.objPath,'String');
saveName = get(handles.saveName,'String');
anonymizeFace   = get(handles.anonymizeFace,'Value');
deleteTxtOutput = get(handles.deleteTxtOutput,'Value');
moveElecInwards = str2double(get(handles.moveElecInwards,'String'));
close

EEG  = evalin('base', 'EEG');
EEG = getChanLocs(EEG, objPath, 'saveName', saveName, 'anonymizeFace', anonymizeFace,...
    'deleteTxtOutput', deleteTxtOutput, 'moveElecInwards', moveElecInwards);
assignin('base','EEG',EEG);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close
