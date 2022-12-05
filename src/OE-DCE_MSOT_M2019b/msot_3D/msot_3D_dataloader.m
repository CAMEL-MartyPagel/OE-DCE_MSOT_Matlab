function varargout = msot_3D_dataloader(varargin)
% MSOT_3D_DATALOADER MATLAB code for msot_3D_dataloader.fig
%      MSOT_3D_DATALOADER, by itself, creates a new MSOT_3D_DATALOADER or raises the existing
%      singleton*.
%
%      H = MSOT_3D_DATALOADER returns the handle to a new MSOT_3D_DATALOADER or the handle to
%      the existing singleton*.
%
%      MSOT_3D_DATALOADER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSOT_3D_DATALOADER.M with the given input arguments.
%
%      MSOT_3D_DATALOADER('Property','Value',...) creates a new MSOT_3D_DATALOADER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msot_3D_dataloader_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msot_3D_dataloader_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msot_3D_dataloader

% Last Modified by GUIDE v2.5 07-Aug-2013 15:35:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msot_3D_dataloader_OpeningFcn, ...
                   'gui_OutputFcn',  @msot_3D_dataloader_OutputFcn, ...
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


% --- Executes just before msot_3D_dataloader is made visible.
function msot_3D_dataloader_OpeningFcn(hObject, eventdata, handles, varargin)

handles.OK = 0;
handles.reconSel = 0;
handles.mspSel = 0;
handles.reconWL = 0;

if (numel(varargin) >= 4 && ~isempty(varargin{1}))
    handles.reconSel = varargin{2};
    handles.reconWL = varargin{3};
    handles.mspSel = varargin{4};
end

if (numel(varargin) >= 1 && ~isempty(varargin{1}))
    handles.datainfo = varargin{1};
    set(handles.e_msotfile,'String',handles.datainfo.XMLFileName);
    refreshLists(handles);
else
    handles.datainfo = [];
end    

guidata(hObject, handles);

uiwait(handles.figure1);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

% --- Outputs from this function are returned to the command line.
function varargout = msot_3D_dataloader_OutputFcn(hObject, eventdata, handles) 
if handles.OK
    varargout{1} = handles.datainfo;
    varargout{2} = handles.reconSel;
    varargout{3} = handles.reconWL;
    varargout{4} = handles.mspSel;
else
    varargout{1} = [];
    varargout{2} = [];
    varargout{3} = [];
    varargout{4} = [];
end
delete(handles.figure1);










function e_msotfile_Callback(hObject, eventdata, handles)
if (~exist(get(hObject,'String'),'file'))
    errordlg('File does not exist');
    return
end
try
    handles.datainfo = loadMSOT(get(hObject,'String'));
    guidata(hObject,handles);
    refreshLists(handles);
catch ex
    errordlg(ex.message);
end

function e_msotfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b_browse_Callback(hObject, eventdata, handles)
[newfile newpath] = uigetfile('*.msot');
if ~isempty(newfile)
    set(handles.e_msotfile,'String',[newpath newfile]);
    e_msotfile_Callback(handles.e_msotfile,[],handles);
end

function refreshLists(handles)
di = handles.datainfo;
if isempty(di), return; end;



rt = cell(numel(di.ReconNode),1);
for j = 1:numel(di.ReconNode)
    rt{j} = di.ReconNode(j).Name;
end
set(handles.lb_recon,'String',rt);
if (handles.reconSel), 
    set(handles.lb_recon,'Value',handles.reconSel); % select last
else
    set(handles.lb_recon,'Value',j); % select last
    handles.reconSel = j;
end
guidata(handles.figure1,handles);
lb_recon_Callback(handles.lb_recon,[],handles);

mt = cell(numel(di.MSPNode),1);
for j = 1:numel(di.MSPNode)
    mt{j} = di.MSPNode(j).Name;
end
set(handles.lb_msp,'String',mt);
if handles.mspSel,
    set(handles.lb_msp,'Value',handles.mspSel); % select last
else
    set(handles.lb_msp,'Value',j); % select last
    handles.mspSel = j;
    guidata(handles.figure1,handles);
end
lb_msp_Callback(handles.lb_msp, [], handles);


function lb_recon_Callback(hObject, eventdata, handles)
ri = get(handles.lb_recon,'Value');
handles.reconSel = ri;
rn = handles.datainfo.ReconNode(ri);
wlarr = unique(cell2mat({handles.datainfo.ScanFrames(rn.ReconStructure(:)).Wavelength}))';

set(handles.l_rname,'String',rn.Name);
set(handles.l_rpixels,'String',num2str(rn.Resolution));
set(handles.l_rroi,'String',num2str(rn.ROI*1e3,'%.1f mm'));
set(handles.s_bgwl,'String',num2str(wlarr,'%i nm'));
if (handles.reconWL)
    set(handles.s_bgwl,'Value',handles.reconWL);
else
    set(handles.s_bgwl,'Value',numel(wlarr));
    handles.reconWL = numel(wlarr);
end
guidata(hObject,handles);


function lb_msp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb_recon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lb_msp_Callback(hObject, eventdata, handles)
mi = get(handles.lb_msp,'Value');
mn = handles.datainfo.MSPNode(mi);
handles.mspSel = mi;
guidata(hObject,handles);

set(handles.l_mmethod,'String',mn.Method);
set(handles.l_mrecon,'String',handles.datainfo.ReconNode(mn.ReconNodeID).Name);
set(handles.l_mwl,'String',num2str(mn.Wavelengths));
set(handles.l_mspec,'String',mn.InputSpectra);













function b_load_Callback(hObject, eventdata, handles)
handles.OK = 1;
guidata(hObject,handles);
close(handles.figure1);

% --- Close without Loading
function b_cancel_Callback(hObject, eventdata, handles)
close(handles.figure1);




function s_bgwl_Callback(hObject, eventdata, handles)
handles.reconWL = get(hObject,'Value');
guidata(hObject,handles);

function s_bgwl_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
