function varargout = msot_diag(varargin)
% MSOT_DIAG MATLAB code for msot_diag.fig
%      MSOT_DIAG, by itself, creates a new MSOT_DIAG or raises the existing
%      singleton*.
%
%      H = MSOT_DIAG returns the handle to a new MSOT_DIAG or the handle to
%      the existing singleton*.
%
%      MSOT_DIAG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSOT_DIAG.M with the given input arguments.
%
%      MSOT_DIAG('Property','Value',...) creates a new MSOT_DIAG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msot_diag_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msot_diag_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msot_diag

% Last Modified by GUIDE v2.5 13-Oct-2014 16:50:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msot_diag_OpeningFcn, ...
                   'gui_OutputFcn',  @msot_diag_OutputFcn, ...
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


% --- Executes just before msot_diag is made visible.
function msot_diag_OpeningFcn(hObject, eventdata, handles, varargin)

% opengl software;

handles.output = hObject;
handles = initialise(handles);
set(handles.f_main,'Name',['msot_diag v1.12']);

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes msot_diag wait for user response (see UIRESUME)
% uiwait(handles.f_main);


% --- Outputs from this function are returned to the command line.
function varargout = msot_diag_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;



















% *************************************************************************
% CREATE FUNCTIONS
% *************************************************************************

function e_indir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lb_setlist_Callback(hObject, eventdata, handles)
function lb_setlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_calib_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lb_tfiles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_geo_serial_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_geo_radius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_geo_coverage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_geo_elements_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_curv_spheres_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_curv_dist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_curv_elements_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb_curv_sim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb_curv_msmts_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% *************************************************************************
% INITIALISATION
% *************************************************************************
function handles = initialise(handles)
    handles.indir = '';
    handles.calibdir = '';
    handles.data = cell([0]);





% *************************************************************************
% INPUT / OUTPUT directories
% *************************************************************************

% --- Browse Input Directory
function b_inbrowse_Callback(hObject, eventdata, handles)
    newpath = uigetdir(handles.indir);
    set(handles.e_indir,'String',newpath);
    load_all(handles);

% --- Edit Input Directory
function e_indir_Callback(hObject, eventdata, handles)
    load_all(handles);

% --- Browse Calibration Directory
function b_calibbrowse_Callback(hObject, eventdata, handles)
    newpath = uigetdir(handles.calibdir);
    set(handles.e_calib,'String',newpath);
    set_calibdir(handles);
    
% --- Edit Calibration Directory    
function e_calib_Callback(hObject, eventdata, handles)
    set_calibdir(handles);








% *************************************************************************
% Data Loading
% *************************************************************************

% --- Load XML Information for Directory
function load_all(handles)
    % check if directory is valid
    indir = get(handles.e_indir,'String');
    if (~exist(indir,'dir'))
        errordlg('Selected Input Folder does not exist');
        return;
    end
    % set directory to handles
    handles.indir = indir;
    % find all XML Files
    cfiles = rdir([indir '\**\*.msot']);
% ** sort based on file date -> not clever
%     [tmp ind] = sort(cell2mat({cfiles.datenum}));
%     cfiles = cfiles(ind);
    % read each file
    d = cell([numel(cfiles),1]);
    for jj = 1:numel(cfiles)
%         fprintf('File %i: %s\n',jj,cfiles(jj).name);
        try
            str = loadMSOT(cfiles(jj).name);
        catch ex
            errordlg(['In: ' cfiles(jj).name '\n' ex.message]);
            return;
        end
        d{jj} = str;
    end
    % sort according to record date
    [tmp ind] = sort(cellfun(@(x) x.CreationTime,d));
    d = d(ind);
    % save to handles array
    handles.data = d;
    guidata(handles.f_main,handles);
    % update GUI Listbox
    update_setlist(handles);
    update_curv_msmts(handles);

 
% --- Check if new sets were added
function b_checknew_Callback(hObject, eventdata, handles)
    % load old and find number
    lastset = handles.data{end};
    maxnum = str2num(strrep(lastset.FolderName,'Scan_',''));
    
    % identify new datasets
    cfiles = rdir([handles.indir '\**\*.msot']);
    add = cell(0,0);
    c=0;
    for j = 1:numel(cfiles)
        cf = cfiles(j).name;
        sep = strfind(cf,'\');
        num = str2num(strrep(cf(sep(end)+6:end),'.msot',''));
        if num > maxnum
            c = c+1;
            add{c} = cfiles(j);
        end
    end
    
    % exit if no new ones found
    if numel(add) == 0,        return;    end
        
    % sort them by date (to add them in the correct order)    
    [tmp ind] = sort(cellfun(@(x) x.datenum,add));
    add = add(ind);
    
    % loadMSOT and add them
    cc = numel(handles.data);
    for jj = 1:numel(add)
        try
            str = loadMSOT(add{jj}.name);
            cc = cc+1;
            handles.data{cc} = str;
        catch ex
            errordlg(['In: ' add{jj}.name '\n' ex.message]);
            return;
        end
    end
    
    % update handles and update list
    guidata(handles.f_main,handles);
    update_setlist(handles);
    update_curv_msmts(handles);
    
    

    
 % --- Set calibration file directory nad update lists
 function handles = set_calibdir(handles)
     % check if directory is valid
    calibdir = get(handles.e_calib,'String');
    if (~exist(calibdir,'dir'))
        errordlg('Selected Calibration Folder does not exist');
        return;
    end
    % set directory to handles
    handles.calibdir = calibdir;
    guidata(handles.f_main,handles);

    update_tcalib(handles);
    update_simfiles(handles);
   
 

% --- Update Listbox with Titles
function handles = update_setlist(handles)
    lb = handles.lb_setlist;
    str = cell([numel(handles.data) 1]);
    % add entry for each file in data
    for jj = 1:numel(handles.data)
        d = handles.data{jj};
        str{jj} = [ d.Name ' (' d.Comment ', ' datestr(d.CreationTime) ')' ];
    end
    set(lb,'String',char(str));



    
    
    
    
% *************************************************************************
% UTILS
% *************************************************************************

function sets = get_selected_sets(handles)
    ind = get(handles.lb_setlist,'Value');
    sets = handles.data(ind);
    











% *************************************************************************
% NOISE
% *************************************************************************


% --- Noise diagram
function b_noisediag_Callback(hObject, eventdata, handles)
sets = get_selected_sets(handles);
for j = 1:numel(sets)
    cset = sets{j};
    par.average = false;
    par.usePower = 0;
    sigMat = loadMSOTSignals(cset,cset.ScanStructure(1,1,1,1,:),par);
    svec = size(sigMat);
    sigMat = double(reshape(sigMat,svec(1),svec(2),prod(svec(3:end))));
    filt = [get(handles.s_highpass,'Value') get(handles.s_lowpass,'Value')]*1e3
%    sigMat(1,:,:,:,:,:) = sigMat(2,:,:,:,:,:);
    sigMat = sigMat(3:end,:,:,:,:,:);
    get_noiselevel(sigMat,50,1,[],[handles.f_main.Name ' - ' cset.Name],filt);
%     get_stddist(sigMat,cset.Name);
end

function s_highpass_Callback(hObject, eventdata, handles)
v = get(hObject,'Value');
set(handles.l_highpass,'String',[num2str(round(v),'%d') ' kHz']);

function s_lowpass_Callback(hObject, eventdata, handles)
v = get(hObject,'Value');
set(handles.l_lowpass,'String',[num2str(round(v),'%d') ' kHz']);











% *************************************************************************
% GEOMETRY
% *************************************************************************

% --- Update Geometry Calibration Files
function update_tcalib(handles)
    tfiles = dir([handles.calibdir '\*.tcal']);
    tlist = cell(1,numel(tfiles));
    for j = 1:numel(tfiles)
        tlist{j} = strrep(tfiles(j).name,'.tcal','');
    end
    set(handles.lb_tfiles,'String',tlist);

    
% --- Geoemtry Selection --> Load tdata file
function lb_tfiles_Callback(hObject, eventdata, handles)
    tfilelist = cellstr(get(hObject,'String'));
    tfile = tfilelist{get(hObject,'Value')};
    load('-mat',[handles.calibdir '\' tfile '.tcal']);
    
    set(handles.e_geo_serial,'String',tfile);
    set(handles.e_geo_radius,'String',num2str(r_sensor*1000,'%.2f'));
    set(handles.e_geo_elements,'String',num2str(numel(angle_sensor),'%i'));
    if exist('coverage','var') == 0
        coverage = uint16(abs(diff(angle_sensor([1 end]))/2/pi*360+diff(angle_sensor([1 2]))/2/pi*360));
    end
    set(handles.e_geo_coverage,'String',num2str(coverage,'%i'));


% --- Find COR
function b_geo_COR_Callback(hObject, eventdata, handles)
% hObject    handle to b_geo_COR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    cset = get_selected_sets(handles);
    if numel(cset) > 1
        errordlg('Only select one scan');
    end
    par = struct();
    par.filter = 2;
    par.f_LPF = 0;
    par.f_HPF = 750e3;
    par.f_HPFOrder = 1;
    par.f_HPFRipple = 1;

    sigMat = loadMSOTSignals(cset{1},[],par);
    sensors = [22 107];
    if strcmp(get(handles.e_geo_elements,'String'),'256')
        sensors = [43 214];
    end
    r_sensor = diag_COR(sigMat,[sensors 1],cset{1}.AverageTemperature);
    set(handles.e_geo_radius,'String',num2str(r_sensor*1e3,'%.2f'));
   
    


% --- Create new tdata file
function b_geo_create_Callback(hObject, eventdata, handles)
    % Coverage
    try coverage = str2num(get(handles.e_geo_coverage,'String'));
    catch ex
        errordlg('Please enter numeric coverage angle'); return;
    end
    
    % Number of elements
    try numel = str2num(get(handles.e_geo_elements,'String'));
    catch ex
        errordlg('Please enter numeric number of elements'); return;
    end
    
    % Radius
    try r_sensor = str2num(get(handles.e_geo_radius,'String'))*1e-3;
    catch ex
        errordlg('Please enter numeric radius'); return;
    end

    % Serial
    serial = get(handles.e_geo_serial,'String');
    if isempty(serial)
        errordlg('Please enter transducer serial number'); return;
    end
    
    % Directory
    if isempty(handles.calibdir) || exist(handles.calibdir,'dir') == 0
        errordlg('Please specify Calibration Directory'); return;
    end
    
    % Check if exists and should be overwritten
    tfilename = [handles.calibdir '\' serial '.tcal']
    if exist(tfilename,'file')
        answer = questdlg('Already exists - Overwrite?');
        if ~strcmp(answer,'Yes')
            return;
        end
    end
    
    imp_resp = [];
    elementstep = coverage/360*2*pi/numel;
    firstelement = (180-coverage)/2/360*2*pi + (coverage/360*2*pi/numel/2);
    angle_sensor = firstelement:elementstep:firstelement+(numel-1)*elementstep;

    save('-mat',tfilename,'coverage','imp_resp','angle_sensor','r_sensor');
    update_tcalib(handles);


function e_geo_serial_Callback(hObject, eventdata, handles)
function e_geo_radius_Callback(hObject, eventdata, handles)
function e_geo_coverage_Callback(hObject, eventdata, handles)
function e_geo_elements_Callback(hObject, eventdata, handles)


    
    
    
    
    
    
% *************************************************************************
% CURVATURE
% *************************************************************************

function e_curv_spheres_Callback(hObject, eventdata, handles)
function e_curv_dist_Callback(hObject, eventdata, handles)
function e_curv_elements_Callback(hObject, eventdata, handles)
function lb_curv_sim_Callback(hObject, eventdata, handles)
function lb_curv_msmts_Callback(hObject, eventdata, handles)




% --- Analyse Curvature
function b_curv_analyse_Callback(hObject, eventdata, handles)
    sets = get_selected_sets(handles);
    if isempty(sets) errordlg('Please select at least one scan'); return; end
        
    % Geometry
    tfilelist = cellstr(get(handles.lb_tfiles,'String'));
    tfilei = get(handles.lb_tfiles,'Value');
    if tfilei == 0 errordlg('Please select geometry'); return; end;
    tfile = [handles.calibdir '\' tfilelist{tfilei} '.tcal']
    
    % Number of spheres
    try spheres = str2num(get(handles.e_curv_spheres,'String'));
    catch ex
        errordlg('Please enter number of spheres'); return;
    end
    
    % Distance
    try dist = str2num(get(handles.e_curv_dist,'String'));
    catch ex
        errordlg('Please enter distance'); return;
    end
    
    % Sensors
    try
        sensors = eval(get(handles.e_curv_elements,'String'));
        sensorsc = get(handles.e_curv_elements,'String');
        sensorsc = strrep(sensorsc,':','-');
    catch ex
        errordlg(['Please enter sensors.' ex.message]);return;
    end
       
    % Simulation
    simfile = [];
    simsel = get(handles.lb_curv_sim,'Value');
    if ~isempty(simsel)
        simlist = get(handles.lb_curv_sim,'String');
        simfile = [handles.calibdir '\' simlist{simsel(1)} '.tsim'];
    end
    
    for jj = 1:numel(sets)
        cset = sets{jj};
        diag_calculate_curvature_BP(cset,spheres,dist,...
            tfile,0,sensors,...
            simfile,[handles.indir '\' cset.Name '_s' sensorsc '.curv'],1);
    end % for loop sets
    
    update_curv_msmts(handles);


 % --- Update Simulation Files
function update_simfiles(handles)
    tfiles = dir([handles.calibdir '\*.tsim']);
    tlist = cell(1,numel(tfiles));
    for j = 1:numel(tfiles)
        tlist{j} = strrep(tfiles(j).name,'.tsim','');
    end
    set(handles.lb_curv_sim,'String',tlist);

    
   
% --- Update Measurement Files
function update_curv_msmts(handles)
    gfiles = dir([handles.indir '\*.curv']);
    if (numel(gfiles)) == 0
        return;
    end
    
    glist = cell(1,numel(gfiles));
    for j = 1:numel(gfiles)
        glist{j} = strrep(gfiles(j).name,'.curv','');
    end
    set(handles.lb_curv_msmts,'String',glist);    



% --- Plot Curvature
function b_curv_plot_Callback(hObject, eventdata, handles)
    sel = get(handles.lb_curv_msmts,'Value');
    if isempty(sel) warndlg('Please select result set'); return; end
    
    glist = cellstr(get(handles.lb_curv_msmts,'String'));
    glist = glist(sel);

    for j = 1:numel(glist)
        glist{j} = [handles.indir '\' glist{j} '.curv'];
    end
    diag_plot_curvature_multiple(glist);




function s_highpass_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function s_lowpass_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
