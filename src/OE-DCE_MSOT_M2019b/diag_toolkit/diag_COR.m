function varargout = diag_COR(varargin)
% DIAG_COR Visualize COR estimation
%
% Parameters:
% - Data description as one of:
%   * filtered (!) sigMat (2D or 3D)
%   * datainfo struct
%   * filename to .msot file as cell (!!!)
% - Array of default values with three items
%   * left channel to use
%   * right channel to use
%   * optional: default slice (3rd dimension in sigMat)

% Last Modified by GUIDE v2.5 13-Jun-2013 09:54:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @diag_COR_OpeningFcn, ...
                   'gui_OutputFcn',  @diag_COR_OutputFcn, ...
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


% --- Executes just before diag_COR is made visible.
function diag_COR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to diag_COR (see VARARGIN)

% Choose default command line output for diag_COR
handles.output = hObject;

if (numel(varargin) == 0)
    errordlg('Insufficient parameters');
    return;
end

% load data
in = varargin{1};
par.filter = true;
par.f_HPF = 200000;          % filter freq

%   if passed as datainfo struct
if (isstruct(in))
    handles.sig = loadMSOTSignals(in,[],par);
    T = in.AverageTemperature;
elseif isnumeric(in)
%   if passed as signal matrix
    handles.sig = in;
    if numel(varargin) >= 3,
        T = varargin{3};
    else
        T = 32.5;
    end
%   if passed as filename
elseif iscell(in)
    [handles.sig di] = loadMSOTSignals(in,[],par);
    T = di.AverageTemperature;
end
handles.c = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;
set(handles.label_c,'String',[num2str(handles.c) ' m/s']);
set(handles.slider_c,'Value',handles.c);

handles.select = [1 2];
handles.slice = 1;
if (numel(varargin) >= 2)
    in2 = varargin{2};
    if (numel(in2) == 3)
        handles.slice = in2(3);
        handles.select = in2(1:2);
    elseif (numel(in2) == 2)
        handles.select = in2;
        handles.slice = 1;
    end
end
set(handles.edit_sel1,'String',num2str(handles.select(1)));
set(handles.edit_sel2,'String',num2str(handles.select(2)));
set(handles.edit_slice,'String',num2str(handles.slice));

handles.COR = 0.0405;
handles.output = handles.COR;
set(handles.label_COR,'String',[num2str(handles.COR) ' m']);
set(handles.slider_COR,'Value',handles.COR);

% Update handles structure
guidata(hObject, handles);
refresh_signal(handles);
uiwait;



% --- Outputs from this function are returned to the command line.
function varargout = diag_COR_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(hObject);


% --- Executes when user attempts to close figure_main.
function figure_main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


% ****** CORE FUNCTIONS *******

function refresh_signal(handles,varargin)
sel = handles.select;
sig = handles.sig(:,:,handles.slice);

r_view = 0.015;
fs = 4e7;
c = handles.c;
COR = handles.COR;
dt = 1/fs;

ax = handles.axes_sig;
if (nargin > 1)
    ax = varargin{1};
end

% first x axis
x1 = (1:size(sig,1))*dt*c;
sig_start = int32( (COR-r_view) * fs/c  ) ;
sig_end = int32( (COR+r_view) * fs/c ) ;
x2 = 2*COR : -dt*c : (2*COR-size(sig,1)*dt*c+dt*c);
sig_start2 = int32( (COR-r_view) * fs/c  ) ;
sig_end2 = int32( (COR+r_view) * fs/c ) ;

axes(ax);
hold off;
plot(x1*1e3,sig(:,sel(1)));
hold all;
% [ymax imax] = find_peaks(sig(sig_start:sig_end,sel(1)),0.95);
% for i = 1:numel(imax)
%     line([imax(i) imax(i)]*dt*c+COR-r_view,ylim,'LineStyle',':','Color','blue');
% end
% if (numel(imax) >= 2)
%     for j = 1:numel(imax)-1
%         text((imax(j)-5)*dt*c+COR-r_view,ymax(j),...
%             ['d_1=' num2str((imax(j+1)-imax(j))*dt*c*1e3,'%.3f') 'mm'],...
%             'Color','blue',...
%             'HorizontalAlignment','left');
% %         ymax
%     end
% end

% if (numel(imax) == 2)
%     text((imax(1)-5)*dt*c+COR-r_view,ymax(1),...
%         ['d_1=' num2str((imax(2)-imax(1))*dt*c*1e3,'%.3f') 'mm'],...
%         'Color','blue',...
%         'HorizontalAlignment','right');
%     ymax
% end

plot(x2*1e3,sig(:,sel(2)));
% plot(x1,sig(:,sel(3)));
% [ymax imax] = find_peaks(sig(sig_start:sig_end,sel(2)),0.25);
% for i = 1:numel(imax)
%     line([x2(sig_start+imax(i)) x2(sig_start+imax(i))],ylim,'LineStyle',':','Color',[0 0.5 0]);
% end
% if (numel(imax) >= 2)
%     for j = 1:numel(imax)-1
%         text(x2(sig_start+imax(j)-5),ymax(j),...
%             ['d_2=' num2str((imax(j+1)-imax(j))*dt*c*1e3,'%.3f') 'mm'],...
%             'Color',[0 0.5 0],...
%             'HorizontalAlignment','left');
% %         ymax
%     end
% end
% 


xlim([COR-r_view COR+r_view]*1e3);
yl = ylim;
line([COR COR]*1e3,ylim,'Color','black','LineStyle',':');
text(COR*1e3+0.2,yl(2)*0.8,...
    sprintf(['COR = ' num2str(COR,'%.5f') ' m\nc = ' num2str(c,'%.1f') ' m/s']),...
    'BackgroundColor','yellow');
hold off;
legend({['Sensor ' num2str(sel(1))],['Sensor ' num2str(sel(2))]},'location','SouthEast');
xlabel('distance from transducer 1 [mm]');
ylabel('OA signal amplitude [a.u.]');
title('Center of Rotation estimation');




% function sigMatFilt = filter_signal(sigMat)
% filter_f = [500 5000] * 1e3;     % lower and higher cut-off frequencies in kHz
% 
% fs = 40e6;                      % 40 Msamples DAQ sampling rate
% 
% f_HPF = filter_f(1) ;
% f_LPF = filter_f(2) ;
% [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .9 ) ;
% [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;
% clear f_HPF f_LPF;
% 
% sigMatFilt = FiltFiltM( b_LPF, a_LPF, sigMat, 1 ) ;
% sigMatFilt = FiltFiltM( b_HPF, a_HPF, sigMatFilt, 1 ) ;
% clear i b_LPF b_HPF a_HPF a_LPF f_HPF f_LPF;



% ****** INTERFACE HANDLERS *******

% --- Executes on slider movement.
function slider_COR_Callback(hObject, eventdata, handles)
set(handles.label_COR,'String',[num2str(get(hObject,'Value')) ' m']);
handles.COR = get(hObject,'Value');
handles.output = handles.COR;
guidata(hObject, handles);
refresh_signal(handles);

% --- Executes on slider movement.
function slider_c_Callback(hObject, eventdata, handles)
set(handles.label_c,'String',[num2str(get(hObject,'Value')) ' m/s']);
handles.c = get(hObject,'Value');
guidata(hObject, handles);
refresh_signal(handles);

function edit_sel1_Callback(hObject, eventdata, handles)
handles.select(1) = str2num(get(hObject,'String'));
guidata(hObject, handles);
refresh_signal(handles);
function edit_sel2_Callback(hObject, eventdata, handles)
handles.select(2) = str2num(get(hObject,'String'));
guidata(hObject, handles);
refresh_signal(handles);
function edit_slice_Callback(hObject, eventdata, handles)
handles.slice = str2num(get(hObject,'String'));
guidata(hObject, handles);
refresh_signal(handles);



% ****** CREATE FUNCTIONS ****** 

function slider_COR_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function slider_c_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function edit_sel1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_sel2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_slice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_newplot.
function button_newplot_Callback(hObject, eventdata, handles)
% hObject    handle to button_newplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
ax = axes();
refresh_signal(handles,ax);


% --- Find Slice
function b_findslice_Callback(hObject, eventdata, handles)
mvec = max(handles.sig(1000:1200,handles.select,:));
figure;plot(squeeze(mvec)');



