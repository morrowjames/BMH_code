function varargout = Fit_Warning(varargin)
% FIT_WARNING MATLAB code for Fit_Warning.fig
%      FIT_WARNING, by itself, creates a new FIT_WARNING or raises the existing
%      singleton*.
%
%      H = FIT_WARNING returns the handle to a new FIT_WARNING or the handle to
%      the existing singleton*.
%
%      FIT_WARNING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIT_WARNING.M with the given input arguments.
%
%      FIT_WARNING('Property','Value',...) creates a new FIT_WARNING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Fit_Warning_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Fit_Warning_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Fit_Warning

% Last Modified by GUIDE v2.5 17-Oct-2016 12:35:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Fit_Warning_OpeningFcn, ...
                   'gui_OutputFcn',  @Fit_Warning_OutputFcn, ...
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


% --- Executes just before Fit_Warning is made visible.
function Fit_Warning_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Fit_Warning (see VARARGIN)

% Choose default command line output for Fit_Warning
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Fit_Warning wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Fit_Warning_OutputFcn(hObject, eventdata, handles) 
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
% CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
delete(handles.figure1);





% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function pushbutton1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
