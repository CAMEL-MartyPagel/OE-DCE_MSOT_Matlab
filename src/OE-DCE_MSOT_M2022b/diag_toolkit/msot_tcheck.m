function varargout = msot_tcheck(varargin)
% MSOT_TCHECK MATLAB code for msot_tcheck.fig

% Last Modified by GUIDE v2.5 19-Jan-2015 10:44:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msot_tcheck_OpeningFcn, ...
                   'gui_OutputFcn',  @msot_tcheck_OutputFcn, ...
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


% --- Executes just before msot_tcheck is made visible.
function msot_tcheck_OpeningFcn(hObject, eventdata, handles, varargin)

warning('off','all');

handles.version = 2;

handles.filename = 'M:\';
handles.datainfo = [];
handles.sigf = [];
handles.recon_cum = [];
handles.peakslice = 0;
handles.c = 0;
handles.r = 0.0405;
handles.roi = 40;
handles.n = 400;
handles.coverage = 270;
handles.numel = 256;
handles.fcenter = 4.9;
handles.bandwidth = 59;
handles.spos = [];
handles.sposa = [];
handles.temp = [];
handles.totsens = [];
handles.sposmap = [];
handles.smod = struct();
handles.pos_tot = 0;
handles.pos_pa = 0;
handles.wbar = [];

handles.impresp = [];
handles.impresp_winsize = 30;
handles.impresp_center = 1000;
handles.impresp_spec = [];
handles.impresp_f = [];
handles.impresp_synth=[];
handles.impresp_synth_spec=[];
handles.impresp_synth_f=[];
handles.impresp_fcenter=[];
handles.impresp_bw=[];
handles.impresp_foc=[];
handles.impresp_foc_spec=[];
handles.impresp_foc_fcenter=[];
handles.impresp_foc_bw=[];

% Simulations from sim file
handles.sim_x = [];
handles.sim_y = [];
handles.sim_Rcurv = [];
handles.sim_sens = [];
handles.sim_fwhm = [];
handles.sim_sens_cum = [];
handles.sim_r = 0.04;
handles.sim_Rcurv_opt = 40;
handles.sim_Rcurv_min = 38;
handles.sim_Rcurv_max = 41;
handles.sim_focus_opt = 38;
handles.sim_focus_min = 36;
handles.sim_focus_max = 39;

% sensor positions
handles.angle_sensor = [];
handles.pos_sensor = [];

% fields from handles that should be saved to file
handles.save_fields = {'smod','filename','datainfo','peakslice','c','r','temp'...
    'spos','sposa','totsens','sposmap','pos_tot','pos_pa','pos_sensor',...
    'impresp','impresp_winsize','impresp_center','impresp_spec','impresp_f',...
    'impresp_synth','impresp_synth_spec','impresp_synth_f','impresp_fcenter','impresp_bw',...
    'impresp_foc','impresp_foc_spec','impresp_foc_fcenter','impresp_foc_bw',...
    'recon_cum','coverage','numel','bandwidth','fcenter','n','roi','angle_sensor',...
    'laserenergy_var','version'};

% Choose default command line output for msot_tcheck
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% *** Load Data - if fails exit
loaddata(handles,[],hObject);

% UIWAIT makes msot_tcheck wait for user response (see UIRESUME)
% uiwait(handles.f_tcheck);


% --- Outputs from this function are returned to the command line.
function varargout = msot_tcheck_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;










% *************************************************************************
% DATA Loading
% *************************************************************************

function ret = loaddata(handles,event,hObject,varargin)
ret = 0;
if numel(varargin) >= 1 && ~isempty(varargin{1}),
    npath = [];
    nfile = varargin{1};
else
    % dialog box, defaults to last file
    [nfile npath] = uigetfile('*.msot','Choose MSOT File',handles.filename);
end

if nfile,
    try
        handles.filename = [npath nfile];
        sep = strfind(nfile,'.');
        nmain = nfile(1:sep(end)-1);
        tcheckfile = [npath nmain '.tcheck'];
        
        if exist(tcheckfile,'file'),
            handles = load_dataset(handles,tcheckfile);
            handles = load_simulations(handles);
            % *** Plot
            plot_cumulative(handles,hObject);
            plot_single_sensor(handles,hObject);
%             set(handles.b_recalc,'Enable','on');
        else
            handles.datainfo = loadMSOT_legacy(handles.filename);
            handles = load_signals(handles);
            handles.temp = handles.datainfo.AverageTemperature;
%             set(handles.b_recalc,'Enable','off');
%             handles = precalc(handles);
%             handles = b_recalc_Callback(hObject, [], handles);
        end
        
        % set GUI
        set(handles.e_browse,'String',handles.filename);
        set(handles.l_name,'String',handles.datainfo.Name);
        set(handles.l_comment,'String',handles.datainfo.Comment);
        set(handles.sl_sensor,'Min',1);
        set(handles.sl_sensor,'Max',size(handles.sigf,2));
        set(handles.sl_sensor,'SliderStep',[1 16]./(size(handles.sigf,2)-1));
        set(handles.e_c,'String',num2str(handles.c));
        set(handles.e_r,'String',num2str(handles.r,'%.5f'));
        set(handles.e_temp,'String',num2str(handles.temp,'%.1f'));
        set(handles.e_peak,'String',num2str(handles.peakslice));
        set(handles.e_elements,'String',num2str(size(handles.sigf,2)));
        set(handles.e_coverage,'String',num2str((handles.datainfo.HWDesc.EndAngle-handles.datainfo.HWDesc.StartAngle)/2/pi*360,'%.1f'));
        
        ret = 1;
    catch ex
        errordlg(['Error loading data: ' ex.message]);
    end
end
guidata(hObject,handles);


% --- Browse Button
function b_browse_Callback(hObject, eventdata, handles)
loaddata(handles,[],hObject);

% --- Enter filename directly
function e_browse_Callback(hObject, eventdata, handles)
loaddata(handles,[],hObject,get(hObject,'String'));


function e_browse_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = load_signals(handles);
% Load Signals (filtered)
par.filter = 2;
par.f_HPF = 200e3;
par.f_HPFOrder = 1;
par.f_HPFRipple = 1;
par.average = 1;
par.useSingle = true;
svec = size(handles.datainfo.ScanStructure);
% if not averaged, then only average the last 25% of shots
% background: energy ramping in DPSS Laser
if (numel(svec) > 4),
    shotstart = round(svec(end)*0.75);
    par.shotsel = shotstart:svec(end);
end
handles.sigf = loadMSOTSignals(handles.datainfo,[],par);

% special fix for SONOTEC miscabling - SM 02.09.2014
if strcmp(handles.datainfo.CreationTimeTxt,'2014-07-30T13:03:24.767+02:00')
    handles.sigf(:,1:32,:,:,:,:) = -handles.sigf(:,1:32,:,:,:,:);
end

handles.pos_tot = handles.datainfo.RunNum;
handles.pos_pa = ceil(sqrt(double(handles.datainfo.RunNum)));

% if dataset incomplete, delete last position to not analyse incomplete
% datapoint
if handles.pos_tot < handles.pos_pa^2,
    handles.pos_tot = handles.pos_tot - 1;
    handles.sigf = handles.sigf(:,:,1:handles.pos_tot,:,:,:);
    warning('Incomplete dataset, truncating last position... Now %i positions',handles.pos_tot);
end



function handles = load_dataset(handles,tcheckfile)

load(tcheckfile,'-mat');
ver = handles.version;
handles = copy_params(handles,saveset);

% we are younger than the version
if ver < handles.version
    warndlg(sprintf('Version of dataset (v%i) newer than application (v%i). Compatibility is not guaranteed, please use current tool version %i.X',handles.version,ver,handles.version));
% we are more mature than the saved dataset
elseif ver > handles.version
    warndlg(sprintf('Version of dataset (v%i) older than application (v%i). Compatibility is not guaranteed, please use current tool version %i.X',handles.version,ver,handles.version));
end
handles.version = ver;
    
handles.sigf = zeros(0,handles.datainfo.HWDesc.NumDetectors);


function dumpAxes(ax,handles,fn,varargin)
f=figure;
set(gcf,'PaperPositionMode','auto');
na = copyobj(ax,f);
if (numel(varargin) >= 1), copyobj(varargin{1},na); end;
for j = 1:numel(na),
    set(na(j),'Units','pixel');
    p = get(na(j),'Position');set(na(j),'Position',[50 50 p(3) p(4)]);
end
set(f,'Position',[100 100 p(3)+100 p(4)+100]);
print(f,'-dpng',[handles.datainfo.StudyPath '\' fn]);
close(f);


















% *************************************************************************
% PRECALCULATIONS
% *************************************************************************

% function b_precalc_Callback(hObject, eventdata, handles)
% precalc(handles);



function handles = precalc(handles)
if isempty(handles.sigf),
    handles = load_signals(handles);
end
    
wbar = waitbar(0,'Precalculations...');

handles = load_simulations(handles);

% Transducer Settings
handles.r = str2num(get(handles.e_r,'String'));
handles = generate_transducer_settings(handles);
handles = generate_synthetic_impresp(handles);

handles = find_peakslice(handles,wbar);
handles = find_sos(handles,wbar);
handles = recon_cum(handles);
handles = plot_cumrecon(handles);

guidata(handles.f_tcheck,handles);
set(handles.b_recalc,'Enable','on');

close(wbar);



function handles = load_simulations(handles);
% Load Simulations
simfile = ['MSOTII_sim' num2str(handles.datainfo.HWDesc.NumDetectors)];
if ~exist([simfile '.mat'],'file'),
    warndlg(sprintf('Simulation File %s does not exist, using default 128el 270° definition',simfile));
    simfile = ['MSOTII_sim128'];
end

load(simfile);
handles.sim_x = x;
handles.sim_y = y;
handles.sim_r = RCOR;
handles.sim_Rcurv = Rcurve_list;
handles.sim_Rcurv_min = Rcurve_min;
handles.sim_Rcurv_opt = Rcurve_opt;
handles.sim_Rcurv_max = Rcurve_max;
handles.sim_sens = sens;
handles.sim_fwhm = fwhm;
handles.sim_sens_cum = sens_cum;
% handles.angle_sensor = angle_sensor(end:-1:1)*2*pi/360;

% Simulated focus
simind = find(handles.sim_Rcurv == handles.sim_Rcurv_opt);
simsens=handles.sim_sens(40,:,simind);
[sx2 ix2] = sort(simsens,'descend');
handles.sim_focus_opt = (mean(handles.sim_x(ix2(1:5)))+handles.sim_r)*1e3;

simind = find(handles.sim_Rcurv == handles.sim_Rcurv_min);
simsens=handles.sim_sens(40,:,simind);
[sx2 ix2] = sort(simsens,'descend');
handles.sim_focus_min = (mean(handles.sim_x(ix2(1:5)))+handles.sim_r)*1e3;

simind = find(handles.sim_Rcurv == handles.sim_Rcurv_max);
simsens=handles.sim_sens(40,:,simind);
[sx2 ix2] = sort(simsens,'descend');
handles.sim_focus_max = (mean(handles.sim_x(ix2(1:5)))+handles.sim_r)*1e3;





%% get peak slice using all sensors and positions
function handles = find_peakslice(handles,wbar)
sigf = handles.sigf;
for pos = 1:size(sigf,3)
    s = squeeze(sigf(500:1600,:,pos,:));
    s = reshape(s,size(s,1)*size(s,2),size(s,3));
    sm = max(s)-min(s);
    [m,i] = max(sm);
    peakslices(pos) = i;
end
handles.peakslice = round(mean(peakslices));
guidata(handles.f_tcheck,handles);
set(handles.e_peak,'String',num2str(handles.peakslice));

%% Transducer Settings
function handles = generate_transducer_settings(handles)
handles.coverage = str2num(get(handles.e_coverage,'String'));
handles.numel = str2num(get(handles.e_elements,'String'));
handles.angle_sensor = linspace(handles.datainfo.HWDesc.EndAngle,handles.datainfo.HWDesc.StartAngle,handles.datainfo.HWDesc.NumDetectors);
handles.pos_sensor = [ cos( handles.angle_sensor' ) sin( handles.angle_sensor' ) zeros( size( handles.angle_sensor' ) ) ] * handles.r ;


%% Focus
function handles = find_sos(handles,wbar);
sigcum = sum(handles.sigf(:,:,:,handles.peakslice),3);
sigcum(1:500,:) = 0;
i = 0;
if (isempty(get(handles.e_temp,'String'))),
    clist = 1460:5:1560;
else
    T = str2double(get(handles.e_temp,'String'));
    clist = round(1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 );
end

rlist = (40.2:0.05:40.8)*1e-3;
% rlist = handles.r;
cc = 0;
waitbar(0,wbar,'Focussing')
for cx = clist
    rc = 0;
    for rx = rlist
        cc = cc + 1;
        waitbar(cc/(numel(rlist)*numel(clist)),wbar);
        rc = rc + 1;
    
        par.c = cx;
        par.r_sensor = rx;
        par.n = 400;
        par.roi = 40e-3;
        par.proj = handles.datainfo.HWDesc.NumDetectors;
        par.filter_f = [500e3 6.5e6];
        par.progress = false;
        par.angle_sensor = handles.angle_sensor;
        par.mirror = 2;
        par.impresp = handles.impresp_synth;
        par.useGPU = 1;
        
        R = -reconMSOT(handles.datainfo,par,1,sigcum);
        val = sort(R(:),'descend');
        metr(cc) = mean(val(val > 4*std(val)));
        rall(cc) = rx;
        call(cc) = cx;
        
%         figure;
%         imagesc(R);
%         axis image;
%         colormap gray;
%         title(num2str(cx));
        
%         fprintf('%i: %06.1f %.5f - %.2f\n',i,cx,rx,p(i));
    
    end
end
[m i] = max(metr);
handles.c = round(call(i));
handles.r = rall(i);
guidata(handles.f_tcheck,handles);
set(handles.e_c,'String',num2str(handles.c,'%i'));
set(handles.e_r,'String',num2str(handles.r,'%.5f'));


% Artificial IR
function handles = generate_synthetic_impresp(handles)
f_center=str2double(get(handles.e_frequency,'String'))*1e6;
fwhm = str2double(get(handles.e_bandwidth,'String'))/100*sqrt(2)*f_center;
phase_shift=0.25*pi;

fs = 40e6; 
w=fwhm/(2*sqrt(2*log(2)));
freq=fs*linspace(0,1,2030);
mag_spec=exp(-((freq-f_center).^2/(2*w^2)));
mag_spec=mag_spec+flipud(mag_spec);
phase_spec=(linspace(0,1,2030)*2*pi)+phase_shift;

iimp=real(ifft(mag_spec.*exp(1i*phase_spec)));
iimp=circshift(iimp',-1015);
% iimp=iimp.*gausswin(2030,100).^4;

% norm_iimp=sqrt(sum(iimp.^2));
iimp=(iimp./max(iimp))';
handles.impresp_synth = iimp;

fsig=fft(iimp);
fsig = abs(fsig)./max(abs(fsig));
handles.impresp_synth_spec = fsig;
handles.impresp_synth_f=[1:2030]*(40/2030)*1e6;




% *** Cumulative Image 
function handles = recon_cum(handles)
pxdist = 3;
sigcum = sum(handles.sigf(:,:,:,handles.peakslice),3);
sigcum(1:500,:) = 0;

clear par;
par.r_sensor = handles.r; % override value in settings to be safe
par.c = handles.c;
par.roi = 40e-3;
par.n = 400;
par.proj = handles.datainfo.HWDesc.NumDetectors;
par.filter_f = [500e3 6.5e6];
par.progress = false;
par.angle_sensor = handles.angle_sensor;
par.mirror = 2;
par.impresp = handles.impresp_synth;
par.useGPU = 1;

handles.recon_cum = -reconMSOT(handles.datainfo,par,1,sigcum);




function e_r_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');


function e_c_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');
function e_peak_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');

function e_c_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_r_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_peak_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_pixel_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');

function e_pixel_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_roi_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');

function e_roi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_coverage_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');
function e_coverage_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_elements_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','off');
function e_elements_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_frequency_Callback(hObject, eventdata, handles)
% set(handles.b_recalc,'Enable','off');
handles = plot_cumulative(handles,hObject);

function e_frequency_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_bandwidth_Callback(hObject, eventdata, handles)
% set(handles.b_recalc,'Enable','off');
handles = plot_cumulative(handles,hObject);


function e_bandwidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function e_temp_Callback(hObject, eventdata, handles)
set(handles.b_recalc,'Enable','on');

function e_temp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







































% *************************************************************************
% CALCULATIONS 
% *************************************************************************

function handles = b_recalc_Callback(hObject, eventdata, handles)

%% Precalculations
handles = precalc(handles);

handles.c = str2num(get(handles.e_c,'String'));
handles.peakslice = str2num(get(handles.e_peak,'String'));
guidata(hObject,handles);

% if isempty(handles.sigf),
%     handles = load_signals(handles);
% end
% 
% 
handles.wbar = waitbar(0,'Recalculating');
handles = find_spheres(handles);
handles.smod = sens_model(handles.datainfo,handles.sigf,handles);
handles = find_impresp(handles);
handles = recon_cum(handles);
handles = plot_laserpower(handles); % to store laser power variations
guidata(hObject,handles);

handles = save_dataset(handles);

% % Waitbar
close(handles.wbar);
handles.wbar = [];

% *** Plot
handles = plot_cumulative(handles,hObject);
handles = print_result(handles);
handles = plot_single_sensor(handles,hObject);









% find sphere positions
function handles = find_spheres(handles);
waitbar(0,handles.wbar,'Finding Grid Positions');
clear par;
par.r_sensor = handles.r; % override value in settings to be safe
par.c = handles.c;
par.roi = 40e-3;
par.n = 400;
par.proj = handles.datainfo.HWDesc.NumDetectors;
par.filter_f = [500e3 6.5e6];
par.progress = false;
par.angle_sensor = handles.angle_sensor;
par.mirror = 2;
par.useGPU = 1;
par.impresp = handles.impresp_synth;

dimvec = linspace(-par.roi/2,par.roi/2,par.n);
handles.spos = nan(handles.pos_pa^2,2);
handles.sposa = nan(handles.pos_pa^2,2);
handles.totsens = [];
handles.sposmap = false(par.n);
% all positions
% figure;
Rall = -reconMSOT(handles.datainfo,par,1,handles.sigf(:,:,:,handles.peakslice));
for pos = 1:size(handles.sigf,3);
    waitbar(pos/size(handles.sigf,3),handles.wbar);
    R = Rall(:,:,pos);
%     imagesc(R);
%     pause(0.25);
    [~,i] = max(R(:));
    
    n = size(R,1);
    handles.spos(pos,2) = mod(i,n);
    handles.spos(pos,1) = ceil(i/n);
    handles.sposa(pos,1) = dimvec(handles.spos(pos,1));
    handles.sposa(pos,2) = dimvec(handles.spos(pos,2));
%     loc = false(size(R));
%     loc(handles.spos(pos,2)-pxdist:handles.spos(pos,2)+pxdist,handles.spos(pos,1)-pxdist:handles.spos(pos,1)+pxdist) = true;
%     Rlim = R(loc);
%     handles.totsens(pos) = max(Rlim(:));
%     handles.sposmap = handles.sposmap | loc;
end
guidata(handles.f_tcheck,handles);









% *** Plot laser energy
function handles = plot_laserpower(handles)
par.roi = 40e-3;
par.n = 400;
waitbar(0,handles.wbar,'Plotting Laser Power');

lp = handles.datainfo.LaserEnergy;
if numel(size(lp)) > 2,
    nshots = size(lp,5);
    navg = round(0.25*nshots);
    lp = mean(lp(:,:,1,1,end-navg:end),5);
end
lpo = max(lp');
handles.laserenergy_var = std(lpo)./mean(lpo)*100;  % variation for report

if isfield(handles.datainfo,'Temperature'),
    tpm = mean(handles.datainfo.Temperature(:,:,1,1,1)');
else
    tpm = ones(numel(lpo),1).*handles.datainfo.AverageTemperature;
end
handles.temp_var = std(tpm) ./ mean(tpm) *100;

axes(handles.ax_laserpower);
[a, h1, h2] = plotyy(1:numel(lpo),lpo,1:numel(lpo),tpm);
set(a(1),'FontSize',7);
set(a(2),'FontSize',7);
ylabel(a(1),'Laser Energy (mJ)');
ylabel(a(2),'Water Temp (°C)');
ymax = double(max(lpo(:))*1.1);
if ymax > 0,
    ylim(a(1),[0 ymax]);
end
out = 'OK';
if handles.laserenergy_var > 5, out = 'WARN'; end;
if handles.laserenergy_var > 10, out = 'FAIL'; end;
text(1,1,sprintf('Laser Variation: %.1f%% - %s',handles.laserenergy_var,out),'Color',h1.Color,'FontSize',8,'VerticalAlignment','bottom','FontWeight','bold');

ylim(a(2),[min(tpm)*0.9 max(tpm)*1.1]);
out = 'OK';
if handles.temp_var > 1, out = 'WARN'; end;
if handles.temp_var > 3, out = 'FAIL'; end;
text(1,ymax,sprintf('Temp Variation: %.1f%% - %s',handles.temp_var,out),'Color',h2.Color,'FontSize',8,'VerticalAlignment','top','FontWeight','bold');


% interpolate laser power on grid
lpo = reshape(lpo,handles.pos_pa,handles.pos_pa);
sax = reshape(handles.sposa(:,1),handles.pos_pa,handles.pos_pa);
say = reshape(handles.sposa(:,2),handles.pos_pa,handles.pos_pa);
[X, Y] = meshgrid(nanmean(sax),nanmean(say'));
dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
[Xi, Yi] = meshgrid(dimvec2,dimvec2);
lpo_all = interp2(X,Y,lpo,Xi,Yi,'spline',0);
axes(handles.ax_laserpower_img);
imagesc(lpo_all);
axis image;
axis off;
dumpAxes(handles.p_laser,handles,'LaserEnergy.png');





%% *** Find impulse response
function handles = find_impresp(handles)
srange = 500:1700;
waitbar(0,handles.wbar,'Finding Average Impulse Response');

handles.impresp_center = 1015; % location of the impulse response
handles.impresp_winsize = 30;    % half window size
npos = 1;%numel(handles.smod.pos_focus{1});
nsens = size(handles.sigf,2);
cc = 0; cct = nsens*npos;
sc = zeros(size(handles.sigf,1),cct);
for sens = 1:nsens
    for pos = handles.smod.pos_focus{sens}(1)
       s = handles.sigf(:,sens,pos,round(handles.smod.slice_mean(sens)));
       cc = cc+1;
       waitbar(cc/cct,handles.wbar);
       [ms, spos] = max(s(srange));
       sc(handles.impresp_center-handles.impresp_winsize:handles.impresp_center+handles.impresp_winsize,cc) = ....
            s(srange(1)+spos-handles.impresp_winsize:srange(1)+spos+handles.impresp_winsize) ...
            .* gausswin(handles.impresp_winsize*2+1,4.0) ;%./ ms;
%         [Y, f, NFFT] = calc_fft(sc(:,cc));
%         scf(:,cc) = abs(Y(1:(NFFT/2+1)));
    end
end
% scx = scf(:,1:npos:end);
% [m i] = max(scx,[],1);
% figure;bar(f(i));

impresp_foc = mean(sc(:,1:end),2);
impresp_foc = impresp_foc ./ max(abs(impresp_foc));
handles.impresp_foc = impresp_foc;

[Y, handles.impresp_f, NFFT] = calc_fft(impresp_foc);
handles.impresp_foc_spec = abs(Y(1:(NFFT/2+1)));

[p i] = max(handles.impresp_foc_spec);
handles.impresp_foc_fcenter = handles.impresp_f(i);
i = find(handles.impresp_foc_spec/p > 0.5);
handles.impresp_foc_bw = diff(handles.impresp_f([i(1) i(end)]))/handles.impresp_foc_fcenter*100/sqrt(2);



% find central position
distfromcenter = sqrt(sum(handles.sposa.^2,2));
[cdist cpos] = min(distfromcenter);

% find peak slice 
clear p2p;
for sl = 1:size(handles.sigf,4),
    p2p(sl) = mean(max(handles.sigf(srange,:,cpos,sl)) - min(handles.sigf(srange,:,cpos,sl)));
end
[pint pslice] = max(p2p);

s = handles.sigf(:,:,cpos,pslice);
sc = zeros(size(s));
handles.impresp_center = 1015; % location of the impulse response
handles.impresp_winsize = 30;    % half window size
for sens = 1:size(s,2)
    waitbar(sens/size(s,2),handles.wbar);
    [mi spos] = max(s(srange,sens));
    sc(handles.impresp_center-handles.impresp_winsize:handles.impresp_center+handles.impresp_winsize,sens) = ....
        s(srange(1)+spos-handles.impresp_winsize:srange(1)+spos+handles.impresp_winsize,sens) ...
        .* gausswin(handles.impresp_winsize*2+1,4.0);
end
impresp = mean(sc');
impresp = impresp ./ max(abs(impresp));
handles.impresp = impresp;

[Y handles.impresp_f NFFT] = calc_fft(impresp);
handles.impresp_spec = abs(Y(1:(NFFT/2+1)));

[p i] = max(handles.impresp_spec);
handles.impresp_fcenter = handles.impresp_f(i);
i = find(handles.impresp_spec/p > 0.5);
handles.impresp_bw = diff(handles.impresp_f([i(1) i(end)]))/handles.impresp_fcenter*100/sqrt(2);







function handles = plot_impresp(handles)
axes(handles.ax_cum_impresp);
hold off;
plot(handles.impresp);
xlim([handles.impresp_center-handles.impresp_winsize*2 handles.impresp_center+handles.impresp_winsize*2]); hold all;
set(gca,'FontSize',8);
% dectest = deconvwnr(s(:,1),impresp',20);
% testsig = s(:,102);
% dectest = deconvwnr(impresp',impresp',0);

% subplot(2,2,2);
% plot(testsig);hold all;
% plot(dectest);
% xlim([cent-csz*5 cent+csz*5]);

axes(handles.ax_cum_impspec);
hold off;
h1 = plot(handles.impresp_f*1e-6,handles.impresp_spec./max(handles.impresp_spec)); 
xlim([0 12]);
set(gca,'FontSize',8);
xlabel('Frequency [MHz]');

axes(handles.ax_cum_impspec);
hold on;
h2 = plot(handles.impresp_synth_f*1e-6,handles.impresp_synth_spec);
h3 = plot(handles.impresp_f*1e-6,handles.impresp_foc_spec./max(handles.impresp_foc_spec));
legend({'center','spec','focus'});
text(handles.impresp_fcenter*1e-6,0.5,sprintf('Fc: %.2fMHz\nBW: %.1f%%',handles.impresp_fcenter*1e-6,handles.impresp_bw),'Color',h1.Color,'FontSize',8,'HorizontalAlign','center');
text(handles.impresp_fcenter*1e-6,0.20,sprintf('Fc: %.2fMHz\nBW: %.1f%%',handles.impresp_foc_fcenter*1e-6,handles.impresp_foc_bw),'Color',h3.Color,'FontSize',8,'HorizontalAlign','center');

axes(handles.ax_cum_impresp);
hold on;
plot(handles.impresp_synth);
plot(handles.impresp_foc);
legend({'center','spec','focus'});

dumpAxes(handles.ax_cum_impspec,handles,'ImpulseResponse.png');















%% *** Calculate sensitivity model
function smod = sens_model(datainfo,sigf,handles)

zspac = mean(diff(datainfo.ZPositions(2:end-1)));
znum = numel(datainfo.ZPositions);
zax = 0:zspac:zspac*(znum-1);
waitbar(0,handles.wbar,'Calculating Sensitivity Field');

srange = 500:1500;
sensor = 106;

f = [];
smod.fwhm = nan(size(sigf,2),handles.pos_pa^2);
smod.sensitivity = zeros(size(sigf,2),handles.pos_pa^2);
smod.err = zeros(size(sigf,2),handles.pos_pa^2);
smod.slice = zeros(size(sigf,2),handles.pos_pa^2);
smod.peakdist = zeros(size(sigf,2),handles.pos_pa^2);
smod.slice_std = zeros(size(sigf,2),1);
smod.slice_mean = zeros(size(sigf,2),1);

% f = figure;
cc = 0;
for sensor = 1:size(sigf,2)
    for pos = 1:size(sigf,3)
        cc = cc + 1;
        waitbar(cc/(size(handles.sigf,3)*size(handles.sigf,2)),handles.wbar);

        sigx = squeeze(sigf(:,sensor,pos,:));
%         sigx = deconvwnr(sigx',handles.impresp,0);
        sigx = sigx';
        sigx = sigx(:,srange)';

%         % ** Peak - 2 - Peak
%         [peak ppeak] = max(sigx);
%         sp2p = peak - min(sigx);
        [peak, ppeak] = max(sigx);
        sp2p = peak;
        

        % reduce noise floor
        noisefloor = mean(sp2p(sp2p < max(sp2p)*0.25));
        if isnan(noisefloor), noisefloor = mean(sp2p(sp2p < max(sp2p)*0.5)); end
        sp2p_ = sp2p - noisefloor;
        sp2p_ (sp2p_ < 0) = 0;

        % *** simplistic
%         [sens zsl] = max(sp2p);
%         mu = zsl * zspac;
%         sigma = nnz(sp2p >= 0.5*sens)*zspac/2.35;
        
        % *** gaussian fit
        [sigma, mu, sens] = mygaussfit(zax,sp2p_,0.3);
        sens = sens + noisefloor;
        zsl = round(mu/zspac);   %peak slice
        if (zsl <= 0 || zsl > numel(zax) || isnan(zsl)),
            sigma = nan;
            mu = nan;
            sens = 0;
            zsl = nan;
            peakpos = nan;
        else            
            peakpos = ppeak(zsl)+ srange(1);
        end

        % *** determin error
        gfun = @(x,mu,sigma) exp(-1 * (x - mu).^2 / ( 2 * sigma.^2 ) );
        mp2p = gfun(zax,mu,sigma);      % get modelled curve
        sp2pn = sp2p ./ sens;           % normalise sensitivity
        if (pos == 41 && sensor == 15) || (sensor == 47 && pos == 41)|| (sensor == 109 && pos == 41), 
            sr2 = round(ppeak(zsl));
            sr2 = sr2-50:sr2+100;
            figure('Position',[100 100 1000 300]);
            subplot(1,2,1);
%             imagesc(zax-zax(zsl),(sr2+srange(1))*2.5e-8*1e6,sigx(sr2,:)); 
%             xlabel('Elevation (mm)'); ylabel('Time (µs)'); colormap jet;
%             title(num2str(sensor,'Sensor: %i'));
%             xlim([-3 3]);
            plot((sr2+srange(1))*2.5e-8*1e6,sigx(sr2,zsl));
            xlabel('Time (µs)');
            title(num2str(sensor,'Sensor: %i'));
            subplot(1,2,2);
            plot(zax-zax(zsl),sp2pn);
            hold all;
            plot(zax-zax(zsl),mp2p); 
            xlim([-3 3]);
            xlabel('Elevation (mm)'); legend({'measured','fitted'},'Location','South');
            title(num2str(pos,'Position: %i')); 
        end
        
        merr = mp2p - sp2pn;            % model error
        merrtot = sqrt(sum(merr.^2))./numel(merr);   % overall modelling error 
%         merrtot2 = sqrt(sum((merr.*abs(zax-mu)).^2))./numel(merr);   % overall modelling error 

        % try to correct bad measurements 
        if isnan(sens) || (merrtot > 0.030 && merrtot/sigma > 0.030)
%             figure;plot(sp2pn); hold all; plot(mp2p);bar(merr);

            % alternative: approximate gaussian...
            [val ind] = sort(sp2p,'descend');
            sens2 = mean([val(1) val(1) val(1) val(2) val(3)]);
            zsl2 = mean([ind(1) ind(1) ind(1) ind(2) ind(3)]);
            mu2 = zsl2*zspac;
            sigma2 = nnz(val > 0.5*sens2)*zspac/2.35*1.2;
            
            mp2p2 = gfun(zax,mu2,sigma2);      % get modelled curve
            sp2pn2 = sp2p ./ sens2;           % normalise sensitivity

            merr2 = mp2p2 - sp2pn2;            % model error
            merrtot2 = sqrt(sum(merr2.^2))./numel(merr2);   % overall modelling error 
            
%             plot(sp2pn2); plot(mp2p2);bar(merr2,'g');
            
            if isnan(sigma) || (merrtot > merrtot2),
                sens = sens2;
                mu = mu2;
                sigma = sigma2;
                merrtot = merrtot2;
            end;
            
        end

        smod.fwhm(sensor,pos) = sigma*2.35;
        smod.sensitivity(sensor,pos) = sens;
        smod.err(sensor,pos) = merrtot;
        smod.slice(sensor,pos) = mu;
        smod.peakdist(sensor,pos) = peakpos;
    end
    
    
    
    % rotate sphere positions to find spheres in FOW
    a = -handles.angle_sensor(sensor);
    R = [cos(a) -sin(a);sin(a) cos(a)];
    sppos = R*handles.sposa';
    pos = find(abs(sppos(2,:)) < 5e-3);
%     posp = find(abs(sppos(2,:)) < 1e-3);
    sdist = sqrt(sum((handles.sposa(pos,:) - ones(numel(pos),1)*handles.pos_sensor(sensor,1:2))'.^2));
    ssens = smod.sensitivity(sensor,pos);

    [sx ix] = sort(ssens,'descend');
    ssens = ssens ./ mean(sx(1:5));
%     smod.focus(sensor) = mean(sdist(ix(1:3)));

    for rc = 1:numel(handles.sim_Rcurv)
        
        % find overlay metric
        handles.smod = smod;
        sens = find_sensfield(handles,sensor);
        smod.fielderr(sensor,rc) = find_overlay(handles,sens,sensor,handles.sim_Rcurv(rc));;

        % restricted number of spheres
        simsens=handles.sim_sens(40,:,rc);
        simsens = simsens ./ max(simsens(:));

        simval = spline(handles.sim_x+handles.r,simsens,sdist);
        smod.serr(sensor,rc) = sqrt(sum((simval - ssens).^2));
        
    end
    
%     [errx erri] = min(smod.fielderr(sensor,:));
    [errx erri] = min(smod.serr(sensor,:));
    smod.focus(sensor) = handles.sim_Rcurv(erri);

%     % sensitivity and fwhm in theoretical focus
    focus_offset = handles.r-handles.sim_Rcurv_opt*1e-3;
    pos_focus = find(abs(sppos(2,:)) < 5e-3 & (sppos(1,:)-focus_offset) <= 5e-3 & (sppos(1,:)-focus_offset) >= -5e-3);
%     [ssv ssi] = sort(smod.sensitivity(sensor,:),'descend');
%     pos_focus = ssi(1:7);
    
    smod.fwhm_focus(sensor) = mean(smod.fwhm(sensor,pos_focus));
%     [sensor, pos_focus;nan, smod.fwhm(sensor,pos_focus)]
    smod.sens_focus(sensor) = mean(smod.sensitivity(sensor,pos_focus));
    smod.pos_focus(sensor) = {pos_focus};

    % estimate view angle (only in focus)
    smu = smod.slice(sensor,pos_focus);
    smu = smu(~isnan(smu));
    smod.slice_mean(sensor) = mean(smu);
    smod.slice_std(sensor) = std(smu);
  
end

sens = smod.sens_focus ./ mean(smod.sens_focus(:));
sens = 20*log10(sens);
brok = (sens < -5);

sf = smod.focus;
sf(brok) = nan;
smod.focus_int = full(smooth2a(sf,round(0.05*size(sigf,2))));

si = smod.slice_mean;
si(brok) = nan;
smod.slice_int = full(smooth2a(si,round(0.05*size(sigf,2))));

% actually not used
smod.sens_focus_int = full(smooth2a(smod.sens_focus,round(0.05*size(sigf,2))));







function handles = save_dataset(handles)

waitbar(0,handles.wbar,'Saving Dataset...');

saveset = struct();
for j = 1:numel(handles.save_fields)
   f = handles.save_fields{j};
   
   saveset = setfield(saveset,f,getfield(handles,f));
end

ind = strfind(handles.filename,'\');
dirname = handles.filename(1:ind(end));
xmlname = handles.filename(ind(end)+1:end);
ind = strfind(xmlname,'.');
setname = xmlname(1:ind-1);
filename = [dirname setname '.tcheck'];

save(filename,'saveset','-v7.3');



















% *************************************************************************
% CUMULATIVE PLOTS
% *************************************************************************

function handles = plot_cumulative(handles,hObject);

handles.wbar = waitbar(0,'Plotting Cumulative');

handles = plot_cumrecon(handles);
handles = plot_cumsens(handles);
handles = plot_cumfwhm(handles);
% handles = plot_cumoverlay(handles);
handles = plot_cumcenter(handles);
handles = plot_cumspecs(handles);
handles = plot_cumslices(handles);
handles = plot_laserpower(handles);
handles = plot_impresp(handles);


handles = print_result(handles);



close(handles.wbar);
handles.wbar = [];

dumpAxes(handles.p_cum,handles,'Summary.png');
datacursormode on;



function handles = plot_cumrecon(handles)
par.roi = 40e-3;
par.n = 400;

R = handles.recon_cum;
axes(handles.ax_cum_image);
hold off;
cax = [0 max(R(:))*0.5];
R_c = imthres(R,cax,hot);
image(R_c); 
% hold on;

set(handles.ax_cum_image,'FontSize',7);
cent = round(par.n / 2);
tickdst = 0.5e-2;
pertick = par.n / par.roi  * tickdst;
ntick = floor(cent / pertick);
tick = cent-ntick*pertick:pertick:cent+ntick*pertick;
set(gca,'XTick',tick);
set(gca,'XTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
set(gca,'YTick',tick);
set(gca,'YTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
grid on;

dumpAxes(handles.ax_cum_image,handles,'ReconstructedImage.png');


%% *** plot cumulative sensitivity
function handles = plot_cumsens(handles)
par.roi = 40e-3;
par.n = 400;
waitbar(0,handles.wbar,'Plotting Sensitivity Field');


% sens = sum(handles.smod.sensitivity);
% sens = reshape(sens,10,10);
% 
% sax = reshape(handles.sposa(:,1),10,10);
% say = reshape(handles.sposa(:,2),10,10);
% [X Y] = meshgrid(mean(sax),mean(say'));
% dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
% [Xi Yi] = meshgrid(dimvec2,dimvec2);
% sens_all = interp2(X,Y,sens,Xi,Yi,'spline',0);
% sens_all = sens_all ./ max(sens_all(:));
% 
% % load sim and find contours
% ind = find(handles.sim_Rcurv == 40);
% sens = handles.sim_sens_cum(:,:,ind);
% [X Y] = meshgrid(handles.sim_x,handles.sim_x);
% dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
% [Xi Yi] = meshgrid(dimvec2,dimvec2);
% sens = interp2(X,Y,sens,Xi,Yi,'spline',0);
% c2 = bwperim(sens > max(sens(:))*0.75);
% c3 = bwperim(sens > max(sens(:))*0.5);
% c6 = bwperim(sens > max(sens(:))*0.25);
sens_all = find_sensfield(handles,0);

% find overlay
[overlay c2 c3 c6] = find_overlay(handles,sens_all,0);
sens_all = round_image(sens_all,14);
handles.sens_overlay = overlay;

% plot
axes(handles.ax_cum_sens);
image(ind2rgb(round(sens_all*64),jet)); hold all;
h = image(ind2rgb(round(c2*64),gray));
set(h,'AlphaData',c2);
h = image(ind2rgb(round(c3*64),gray));
set(h,'AlphaData',c3);
h = image(ind2rgb(round(c6*64),gray));
set(h,'AlphaData',c6);
text(5,15,sprintf('Total Field Deviation: %.2f',overlay*100),'Color','white','FontWeight','bold','FontSize',8);

axis on;
set(gca,'FontSize',7);
cent = round(par.n / 2);
tickdst = 0.5e-2;
pertick = par.n / par.roi  * tickdst;
ntick = floor(cent / pertick);
tick = cent-ntick*pertick:pertick:cent+ntick*pertick;
set(gca,'XTick',tick);
set(gca,'XTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
set(gca,'YTick',tick);
set(gca,'YTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
grid on;

dumpAxes(handles.ax_cum_sens,handles,'SensitivityMap.png');







%% *** plot cumulative FWHM
function handles = plot_cumfwhm(handles)
par.roi = 40e-3;
par.n = 400;
waitbar(0,handles.wbar,'Plotting FWHM Field');

% get mean without nans
fwhm = abs(handles.smod.fwhm);
nloc = ~isnan(fwhm); fwhm(~nloc) = 0;
fwhm = sum(fwhm) ./ sum(nloc);
fwhm(isnan(fwhm)) = 0;

% fwhm_std = std(abs(handles.smod.fwhm));
nax = sqrt(size(handles.sposa,1));
fwhm = reshape(fwhm,nax,nax);
% fwhm_std = reshape(fwhm_std,10,10);

sax = reshape(handles.sposa(:,1),nax,nax);
say = reshape(handles.sposa(:,2),nax,nax);
[X Y] = meshgrid(nanmean(sax),nanmean(say'));
dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
[Xi Yi] = meshgrid(dimvec2,dimvec2);
fwhm_all = interp2(X,Y,fwhm,Xi,Yi,'spline',0);
fwhm_all = round_image(fwhm_all,14,[0.9 3]);

axes(handles.ax_cum_fwhm);
fwhm_t = (fwhm_all - 0.9) ./ 2.1;
fwhm_t(fwhm_t < 0) = 0;
fwhm_c = ind2rgb(round(fwhm_t*64),jet);
image(fwhm_c);
axis off;

axis on;
set(gca,'FontSize',7);
cent = round(par.n / 2);
tickdst = 0.5e-2;
pertick = par.n / par.roi  * tickdst;
ntick = floor(cent / pertick);
tick = cent-ntick*pertick:pertick:cent+ntick*pertick;
set(gca,'XTick',tick);
set(gca,'XTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
set(gca,'YTick',tick);
set(gca,'YTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
grid on;

par.cbarfontsize = 7;
par.cbarlabel={'0.9mm','3.0mm'};
plotColorbar(handles.ax_cum_fwhm_cb,jet,par)

dumpAxes(handles.ax_cum_fwhm,handles,'ElevationalResolutionMap.png');











function handles = plot_cumspecs(handles)
waitbar(0,handles.wbar,'Evaluating Specs');

nel = numel(handles.angle_sensor);
nstep = 32;
if nel > 256, nstep = 64; end;

% for sensor = 1:nel
%     fwe = handles.smod.fwhmerr(sensor,:);
%     [merr mind] = min(fwe);
%     simsel(sensor) = mind;
%     
%     sens_all = find_sensfield(handles,sensor);
%     overlay(sensor) = find_overlay(handles,sens_all,sensor);
% %     p = polyfit(handles.sim_Rcurv,log(handles.smod.fwhmerr(sensor,:)),1);
% %     
% %     figure;plot(handles.sim_Rcurv,handles.smod.fwhmerr(sensor,:));
% %     hold all;
% %     plot(handles.sim_Rcurv,polyval(p,handles.sim_Rcurv));
% end
for sensor = 1:nel,
    fwe = handles.smod.serr(sensor,:);
    [merr mind] = nanmin(fwe);
    simsel(sensor) = mind;
end
simind = find(handles.sim_Rcurv == handles.sim_Rcurv_opt);

sens = handles.smod.sens_focus;
msens = nanmean(sens);
sens = sens ./ msens;
sens = 20*log10(sens);
brok = (sens < -5);
xf = 1:nel;
sf = handles.smod.focus_int;
sf(brok) = nan;

axes(handles.ax_cum_specs);
hold off;
% [a h1 h2] = plotyy(1:nel,ones(1,nel)*1.5,1:nel,handles.sim_Rcurv(simsel));
[a h1 h2] = plotyy(1:nel,ones(1,nel)*1.5,xf,sf);
set(a,'FontSize',7);
title('Main Beam');
xlabel('Sensor Elements');
ylabel('Field Deviation');
ylabel(a(2),'Focal Spot (mm)');
set(a,'XLim',[0 nel]);
set(a,'XTick',[1 nstep:nstep:nel]);

% FWHM plot
set(h1,'LineStyle',':');
% set(a(1),'XLim',[0 128]);
set(a(1),'YLim',[0 2]);
set(a(1),'YTick',[0:0.5:4]);

% set(a(2),'XLim',[0 128]);
set(a(2),'YLim',[min(handles.sim_Rcurv)-1 max(handles.sim_Rcurv)+1]);
hold all;
% plot(1:nel,handles.smod.focus*1e3);
h3 = bar(handles.smod.serr(:,simind),'EdgeColor','none','FaceColor',get(h1,'Color'));
h4 = bar(handles.smod.serr(:,simind).*double(brok'),'EdgeColor','none','FaceColor',[0 0 0 ]./255);
% h3 = bar(handles.smod.fielderr(:,simind)*100,'EdgeColor','none','FaceColor',get(h1,'Color'));
line([0 nel],[1 1]*handles.sim_Rcurv_opt,'Color',get(h2,'Color'),'Parent',a(2));
line([0 nel],[1 1]*handles.sim_Rcurv_max,'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
line([0 nel],[1 1]*handles.sim_Rcurv_min,'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
line([0 nel],[1 1]*mean(handles.smod.focus),'Color',get(h2,'Color'),'Parent',a(2),'LineStyle','--');
line([0 nel],[1 1]*mean(handles.smod.serr(:,simind)),'Color',get(h1,'Color'),'Parent',a(1),'LineStyle','--');
% line([0 nel],[1 1]*mean(handles.smod.fielderr(:,simind)),'Color',get(h1,'Color'),'Parent',a(1),'LineStyle','--');
% line([0 nel],[38 38],'Color',get(h2,'Color'),'Parent',a(2));
% line([0 nel],[39 39],'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
% line([0 nel],[36 36],'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
text(0,0,sprintf(' Average Focus: %.1fmm',mean(handles.smod.focus(~brok))),'Color','white','FontWeight','bold','FontSize',8,'VerticalAlignment','bottom');

misshaped = nnz(round(sf*10)/10 > handles.sim_Rcurv_max) + ...
        nnz(round(sf*10)/10 < handles.sim_Rcurv_min);
misshaped_p = misshaped ./ numel(handles.smod.focus)*100;
text(0,2,sprintf(' Misshaped Elements: %i (%.1f%%)',misshaped,misshaped_p),'Color',h2.Color,'FontWeight','bold','FontSize',8,'VerticalAlignment','top');

dumpAxes(fliplr(a),handles,'ElementFocussing.png');





% 
% function handles = plot_cumoverlay(handles)
% waitbar(0,handles.wbar,'Evaluating Specs');
% 
% nel = numel(handles.angle_sensor);
% nstep = 32;
% if nel > 256, nstep = 64; end;
% 
% % for sensor = 1:nel
% %     fwe = handles.smod.fwhmerr(sensor,:);
% %     [merr mind] = min(fwe);
% %     simsel(sensor) = mind;
% %     
% %     sens_all = find_sensfield(handles,sensor);
% %     overlay(sensor) = find_overlay(handles,sens_all,sensor);
% % %     p = polyfit(handles.sim_Rcurv,log(handles.smod.fwhmerr(sensor,:)),1);
% % %     
% % %     figure;plot(handles.sim_Rcurv,handles.smod.fwhmerr(sensor,:));
% % %     hold all;
% % %     plot(handles.sim_Rcurv,polyval(p,handles.sim_Rcurv));
% % end
% for sensor = 1:nel,
%     fwe = handles.smod.fielderr(sensor,:);
%     [merr mind] = nanmin(fwe);
%     simsel(sensor) = mind;
% end
% 
% axes(handles.ax_cum_extrema);
% hold off;
% [a h1 h2] = plotyy(1:nel,ones(1,nel)*2.5,1:nel,handles.sim_Rcurv(simsel));
% set(a,'FontSize',7);
% title('Total Field');
% xlabel('Sensor Elements');
% ylabel('Field Deviation');
% ylabel(a(2),'Curvature (mm)');
% set(a,'XLim',[0 nel]);
% set(a,'XTick',[1 nstep:nstep:nel]);
% 
% % FWHM plot
% set(h1,'LineStyle',':');
% % set(a(1),'XLim',[0 128]);
% set(a(1),'YLim',[0 7]);
% set(a(1),'YTick',[0:5]);
% 
% % set(a(2),'XLim',[0 128]);
% set(a(2),'YLim',[min(handles.sim_Rcurv)-1 max(handles.sim_Rcurv)+1]);
% hold on;
% h3 = bar(handles.smod.fielderr(:,5)*100,'EdgeColor','none','FaceColor',get(h1,'Color'));
% line([0 nel],[1 1]*handles.sim_Rcurv_opt,'Color',get(h2,'Color'),'Parent',a(2));
% line([0 nel],[1 1]*handles.sim_Rcurv_max,'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
% line([0 nel],[1 1]*handles.sim_Rcurv_min,'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
% line([0 nel],[1 1]*mean(handles.sim_Rcurv(simsel)),'Color',get(h2,'Color'),'Parent',a(2),'LineStyle','--');
% line([0 nel],[1 1]*mean(handles.smod.fielderr(:,5)*100),'Color',get(h1,'Color'),'Parent',a(1),'LineStyle','--');
% 




%% *** plot cumcenter // actually focus
function handles = plot_cumcenter(handles);
waitbar(0,handles.wbar,'Plotting Center');
nel = size(handles.sigf,2);
nstep = 32;
if nel > 256, nstep = 64; end;

fwhm = handles.smod.fwhm_focus;
sens = handles.smod.sens_focus;
msens = nanmean(sens);
sens = sens ./ msens;
sens = 20*log10(sens);
brok = (sens < -5);
sens(brok) = nan;

% PLOT
axes(handles.ax_cum_center);
hold off;
[a h1 h2] = plotyy(1:nel,ones(1,nel)*1.1,1:nel,sens);
set(a,'FontSize',7);
title(sprintf('At Focus (%imm)',handles.sim_Rcurv_opt));
xlabel('Sensor Elements');
ylabel(a(2),'Sensitivity (dBm)');
ylabel(a(1),'Elevational Resolution (mm)');

set(a,'XLim',[0 nel]);
set(a,'XTick',[1 nstep:nstep:nel]);

% FWHM plot
set(h1,'LineStyle',':');
% set(a(1),'XLim',[0 128]);
set(a(1),'YLim',[0 3]);
set(a(1),'YTick',[0 1 2 3]);

% set(a(2),'XLim',[0 128]);
set(a(2),'YLim',[-7 7]);
set(a(2),'YTick',[-7 -5 -3 0 3 5 7]);
line([0 nel],[0 0],'Color',get(h2,'Color'),'Parent',a(2));
line([0 nel],[5 5],'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
line([0 nel],[-5 -5],'Color',get(h2,'Color'),'Parent',a(2),'LineStyle',':');
line([0 nel],[1 1]*mean(fwhm(~brok)),'Color',get(h1,'Color'),'Parent',a(1),'LineStyle','--');
hold on;

h3 = bar(fwhm,'EdgeColor','none','FaceColor',get(h1,'Color'));
h4 = bar(fwhm.*double(brok),'EdgeColor','none','FaceColor',[0 0 0]./255);
text(0,0,sprintf(' Mean Elevational Resolution: %.2fmm',mean(fwhm(~brok))),'Color','white','FontWeight','bold','FontSize',8,'VerticalAlignment','bottom');
text(0,3,sprintf(' Maximum Deviation: %.1fdBm\n Broken Elements: %i',nanmax(abs(sens(~brok))),nnz(brok)),'Color',h2.Color,'FontWeight','bold','FontSize',8,'VerticalAlignment','top');

dumpAxes(fliplr(a),handles,'Sensitvity_ElevationalResolution.png');




function handles = plot_cumslices(handles)
waitbar(0,handles.wbar,'Plotting View Angle');
nel = size(handles.sigf,2);
nstep = 32;
if nel > 256, nstep = 64; end;

sens = handles.smod.sens_focus;
msens = nanmean(sens);
sens = sens ./ msens;
sens = 20*log10(sens);
brok = (sens < -5);
sf = handles.smod.slice_int;
sf(brok) = nan;

axes(handles.ax_cum_view);
hold off;
% [a h1 h2] = plotyy(1:nel,ones(1,nel)*0.5,1:nel,handles.smod.slice_mean - mean(handles.smod.slice_mean));
[a h1 h2] = plotyy(1:nel,ones(1,nel)*0.5,1:nel,sf - mean(handles.smod.slice_int(~brok)));
hold on;
set(a,'FontSize',7);
title('Peak Slice (View Angle)');
xlabel('Sensor Elements');
ylabel('STD (mm)');
ylabel(a(2),'Deviation (mm)');
set(a,'XLim',[0 nel]);
set(a,'XTick',[1 nstep:nstep:nel]);

set(h1,'LineStyle',':');
% set(a(1),'XLim',[0 128]);
set(a(1),'YLim',[0 0.4]);
set(a(1),'YTick',[0 0.1 0.2 0.3 0.4]);
set(a(2),'YLim',[-0.5 0.5]);
line([1 nel],[1 1]*0.2,'LineStyle',':','Color',get(h2,'Color'),'Parent',a(2));
line([1 nel],[1 1]*-0.2,'LineStyle',':','Color',get(h2,'Color'),'Parent',a(2));

bar(handles.smod.slice_std,'EdgeColor','none','FaceColor',get(h1,'Color'));
bar(handles.smod.slice_std.*double(brok'),'EdgeColor','none','FaceColor',[0 0 0]./255);
% text(0,0.4,sprintf(' Max View Angle Deviation: %.2fmm',max(abs(handles.smod.slice_mean - mean(handles.smod.slice_mean)))),'Color',h2.Color,'FontWeight','bold','FontSize',8,'VerticalAlignment','top');
text(0,0.4,sprintf(' Max View Angle Deviation: %.2fmm',nanmax(abs(sf - mean(handles.smod.slice_int(~brok))))),'Color',h2.Color,'FontWeight','bold','FontSize',8,'VerticalAlignment','top');

dumpAxes(fliplr(a),handles,'ViewAngle.png');




function sens_all = find_sensfield(handles,sensor);
par.roi = 40e-3;
par.n = 400;

if sensor > 0,
    sens = handles.smod.sensitivity(sensor,:);
else
    sens = nansum(handles.smod.sensitivity);
end
nax = sqrt(size(handles.sposa,1));
sens = reshape(sens,nax,nax);


% interpolate sensitivity field
sax = reshape(handles.sposa(:,1),nax,nax);
say = reshape(handles.sposa(:,2),nax,nax);
[X Y] = meshgrid(nanmean(sax),nanmean(say'));
sens(isnan(sens)) = 0;
dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
[Xi Yi] = meshgrid(dimvec2,dimvec2);
sens_all = interp2(X,Y,sens,Xi,Yi,'spline',0);
% sens_all = sens_all - min(sens_all(:));
sens_all = sens_all ./ max(sens_all(:));




function [overlay c2 c3 c6] = find_overlay(handles,sens_all,sensor,varargin);

curv = handles.sim_Rcurv_opt;
if numel(varargin) >= 1,
    curv = varargin{1};
end;

par.roi = 40e-3;
par.n = 400;

% angles
angle_deg = handles.angle_sensor/2/pi*360;

% load sim and find contours
ind = find(handles.sim_Rcurv == curv);
if sensor > 0,
    sens = imrotate(handles.sim_sens(:,:,ind),180-angle_deg(sensor),'nearest','crop');
else
    sens = handles.sim_sens_cum(:,:,ind);
end

[X Y] = meshgrid(handles.sim_x,handles.sim_x);
dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
[Xi Yi] = meshgrid(dimvec2,dimvec2);
sens = interp2(X,Y,sens,Xi,Yi,'spline',0);
c2 = bwperim(sens > max(sens(:))*0.75);
c3 = bwperim(sens > max(sens(:))*0.5);
c6 = bwperim(sens > max(sens(:))*0.25);

% find overlay
sens(sens_all == 0) = 0;
sens = sens ./ max(sens(:));
sens_all2 = sens_all ./ max(sens_all(:));
maxoverlay = sum(sum(sens.*sens));
% overlay = sum(sum(abs(sens-sens_all2).*sens.*sens));
% overlay = sum(sum(abs(sens_all-sens)./sens))./maxoverlay;
overlay = abs((sens_all-sens).*sens);
% figure;
% subplot(1,3,1);
% imagesc(overlay);
% axis image;
% colorbar;
% subplot(1,3,2);
% imagesc(sens_all);
% axis image;
% colorbar;
% title('Meas');
% subplot(1,3,3);
% imagesc(sens);
% axis image;
% colorbar;
% title('Sim');


loc = isnan(overlay);
overlay = overlay(~loc);
overlay = mean(overlay(:));

















%% *** REPORT

function handles = print_result(handles)
res = {};

% Laser Energy Variations
% out = 'PASS';
% if handles.laserenergy_var > 5, out = 'WARN'; end;
% if handles.laserenergy_var > 10, out = 'FAIL'; end;
% res(1,:) = {'Laser Energy Variation (%)',num2str(handles.laserenergy_var,'%.1f'),out,' ','5','10'};
% res.Properties.VariableNames = {'Criterion','Meas','Spec','Warn','Crit','Acceptance'};

% Radius
out = 'PASS';
if abs(handles.r-0.0405) >= 0.0002, out = 'WARN'; end;
if abs(handles.r-0.0405) >= 0.0005, out = 'FAIL'; end;
res(1,:) = {'Radius (mm)',num2str(handles.r*1e3,'%.2f'),out,'40.50','0.20','0.50'};

% Impulse Reponse
fcenter = handles.impresp_fcenter*1e-6;
fc_spec = str2double(get(handles.e_frequency,'String'));
out = 'PASS';
if abs(fcenter-fc_spec) >= 0.2, out = 'WARN'; end;
if abs(fcenter-fc_spec) >= 0.5, out = 'FAIL'; end;
% res(2,:) = {'Center Frequency (MHz)',num2str(fcenter,'%.2f'),out,num2str(fc_spec,'%.2f'),'0.2','0.5'};
res(2,:) = {'Center Frequency (MHz)',num2str(fcenter,'%.2f'),'',num2str(fc_spec,'%.2f'),'',''};

bw = handles.impresp_bw;
out = 'PASS';
if bw <= 60, out = 'WARN'; end;
if bw <= 55, out = 'FAIL'; end;
% res(3,:) = {'Bandwidth (%)',num2str(bw,'%.1f'),out,'> 55','< 60','< 55'};
res(3,:) = {'Bandwidth (%)',num2str(bw,'%.1f'),'','> 55','',''};

% Max View Angle
va = abs(handles.smod.slice_int - mean(handles.smod.slice_int));
% out = ' ';
out = 'PASS';
if max(va) > 0.15, out = 'WARN'; end;
if max(va) > 0.2, out = 'FAIL'; end;
res(4,:) = {'View Angle (mm)',num2str(max(va),'%.2f'),out,'< 0.2mm','> 0.15','>= 0.20'};

% Sensitivity
sens = handles.smod.sens_focus;
msens = nanmean(sens);
sens = sens ./ msens;
sens = 20*log10(sens);
% broken elements --> <-5
brok = sens < -5;
broken = nnz(brok); broken_p = broken ./ numel(sens) *100;
out = 'PASS';
if broken_p > 1, out = 'WARN'; end;
if broken_p >= 2, out = 'FAIL'; end;
res(5,:) = {'Broken Elements',sprintf('%i (%.1f%%)',broken,broken_p),out,'< 2%','> 1%','>= 2%'};


sens(brok) = 0;         % don't count broken for sensitivity variation
out = 'PASS';
if max(abs(sens))>3, out = 'WARN'; end;
if max(abs(sens))>5, out = 'FAIL'; end;
res(6,:) = {'Sensitivity Variation (max, dBm)',num2str(max(abs(sens)),'%.1f'),out,'< 5','> 3','>= 5'};

% Elev Res (focus)
fwhm = mean(handles.smod.fwhm_focus(~brok));
out = 'PASS';
if fwhm>1.1, out = 'WARN'; end;
if fwhm>1.2, out = 'FAIL'; end;
res(7,:) = {'Elevational Resolution (mean, mm)',num2str(fwhm,'%.2f'),out,'< 1.2','> 1.1','>= 1.2'};

% Focus
mfocus = mean(handles.smod.focus);
% out = ' ';
out = 'PASS';
if mfocus < handles.sim_Rcurv_min, out = 'WARN'; end;
if mfocus > handles.sim_Rcurv_max, out = 'FAIL'; end;
res(8,:) = {'Average Focus (mm)',num2str(mfocus,'%.1f'),out,num2str(handles.sim_Rcurv_opt,'%.1f'),num2str(handles.sim_Rcurv_min,'%.1f'),num2str(handles.sim_Rcurv_max,'%.1f')};

% Nb of Misshaped Elements
misshaped = nnz(round(handles.smod.focus_int*10)/10 > handles.sim_Rcurv_max) + ...
        nnz(round(handles.smod.focus_int*10)/10 < handles.sim_Rcurv_min);
misshaped_p = misshaped ./ numel(handles.smod.focus_int)*100;
out = 'PASS';
if misshaped_p > 2, out = 'WARN'; end;
if misshaped_p >= 5, out = 'FAIL'; end;
res(9,:) = {'Misshaped Elements',sprintf('%i (%.1f%%)',misshaped,misshaped_p),out,'< 5%','> 2%','>= 5%'};

% Total Field deviation
out = ' ';
% out = 'PASS';
% if handles.sens_overlay > 1, out = 'WARN'; end;
% if handles.sens_overlay > 1.5, out = 'FAIL'; end;
res(10,:) = {'Total Field Deviation',num2str(handles.sens_overlay*100,'%.2f'),out,' ',' ',' '};

handles.t_result.Data = res;






















% *************************************************************************
% SINGLE ELEMENT
% *************************************************************************

function handles = plot_single_sensor(handles,hObject)

handles = plot_sing_loc(handles);
handles = plot_sing_focus(handles);
handles = plot_sing_sens(handles);
handles = plot_sing_fwhm(handles);


function handles = plot_sing_loc(handles)
sensor = get(handles.sl_sensor,'Value');


axes(handles.ax_sing_loc);
hold off;
plot(handles.pos_sensor(:,1),-handles.pos_sensor(:,2),'.k');
hold on;
plot(handles.pos_sensor(sensor,1),-handles.pos_sensor(sensor,2),'.r','LineWidth',3);
axis off;


function handles = plot_sing_focus(handles);
sensor = get(handles.sl_sensor,'Value');

% load simulation
if handles.numel == 256,
    simfile = 'imasonic_256_15.tsim';
elseif handles.numel == 64,
    simfile = 'imasonic_64_15.tsim';
else
    simfile = 'imasonic_128_14.tsim';
end
load('-mat',simfile);
zres = 0.1;
rsim = handles.sim_r;

% angles
angle_sensor = handles.angle_sensor;
r_sensor = handles.pos_sensor;

% rotate sphere positions to find spheres in FOW
a = -angle_sensor(sensor);
R = [cos(a) -sin(a);sin(a) cos(a)];
sppos = R*handles.sposa';
pos = find(abs(sppos(2,:)) < 5e-3);
% posp = find(abs(sppos(2,:)) < 1e-3);
spos = r_sensor(sensor,1:2);
sdist = sqrt(sum((ones(numel(pos),1)*spos - handles.sposa(pos,:)).^2,2));



axes(handles.ax_sing_loc);
plot(handles.sposa(:,1),-handles.sposa(:,2),'.k');hold on;
plot(handles.sposa(pos,1),-handles.sposa(pos,2),'.','Color',[0 1 1]);
% for j = 1:numel(posp)
% plot(handles.sposa(posp,1),-handles.sposa(posp,2),'.','Color',[0 0 1]);

% focus_offset = handles.r-handles.sim_Rcurv_opt*1e-3;
% pos_focus = find(abs(sppos(2,:)) < 5e-3 & (sppos(1,:)-focus_offset) <= 1e-3 & (sppos(1,:)-focus_offset) >= -2e-3);
pos_focus = handles.smod.pos_focus{sensor};
plot(handles.sposa(pos_focus,1),-handles.sposa(pos_focus,2),'.','Color',[0 0 1]);
%     sdist(j)
%     pause
% end


% distances and measurements
% sdist = sqrt(sum((sppos(:,pos)+(ones(numel(pos),1)*[handles.r 0])').^2));
fwhm = handles.smod.fwhm(sensor,pos);
fwhm_focus = handles.smod.fwhm(sensor,pos_focus);
sdist_focus = sqrt(sum((ones(numel(pos_focus),1)*spos - handles.sposa(pos_focus,:)).^2,2));

sens = handles.smod.sensitivity(sensor,pos);
[sx ix] = sort(sens,'descend');
sens = sens ./ mean(sx(1:5));
sens_focus = handles.smod.sensitivity(sensor,pos_focus);
sens_focus = sens_focus ./ mean(sx(1:5));

[errx erri] = min(handles.smod.serr(sensor,:));
foc = handles.sim_Rcurv(erri);
% foc = mean(sdist(ix(1:3)));
% p = polyfit(sdist,sens',4);
% sensmod = polyval(p,sdist);
% [sfoc ifoc] = max(sensmod);
% foc = sdist(ifoc);


simind = find(handles.sim_Rcurv == handles.sim_Rcurv_opt);
simsens=handles.sim_sens(40,:,simind);
simsens = simsens ./ max(simsens);
simfoc = handles.sim_Rcurv_opt;

simval = spline(handles.sim_x+handles.r,simsens,sdist');
err = handles.smod.serr(sensor,erri);

% plot FWHM
axes(handles.ax_sing_curv);
hold off;
area((simx+rsim)*1e3,max([fwhm_min*zres; fwhm_max*zres]),'FaceColor','green','EdgeColor','none','DisplayName','Tolerance');
set(gca,'FontSize',7);
hold on;
area((simx+rsim)*1e3,min([fwhm_min*zres; fwhm_max*zres]),'FaceColor','white','EdgeColor','none','DisplayName',' ');
plot((simx+rsim)*1e3,fwhm_opt*zres,'Color',[0 0.5 0],'DisplayName','Specification');
plot(sdist*1e3,fwhm,'.k');
plot(sdist_focus*1e3,fwhm_focus,'xk','markerSize',5);
xlim([30 50]);
ylim([0 4]);
ylabel('Slice Thickness FWHM (mm)');
xlabel('Distance from transducer (mm)');
title('Slice Thickness (main beam)');

% plot Sensitivity field
axes(handles.ax_sing_specs);
hold off;
a1 = plot(sdist*1e3,sens,'.'); hold on;
set(gca,'FontSize',7);
a2 = plot((handles.sim_x+handles.r)*1e3,simsens);
a3 = plot(sdist_focus*1e3,sens_focus,'x','MarkerEdgeColor',a1.Color,'MarkerSize',5);
% plot(sdist*1e3,sensmod);
xlim([30 50]);
ylim([0 1.2]);
ylabel('Sensitivity (normalised)');
xlabel('Distance from transducer (mm)');
line([1 1]*handles.r*1e3,ylim,'Color',[0.7 0.7 0.7],'LineStyle',':');
line([1 1]*simfoc,ylim,'Color',get(a2,'Color'),'LineStyle','--');
line([1 1]*foc,ylim,'Color',get(a1,'Color'),'LineStyle','--');
text(foc+0.5,1.1,[num2str(foc,'%.1f') 'mm (' num2str((foc-simfoc),'%+.1f') 'mm)'],'Color',get(a1,'Color'),'FontSize',7);
legend({['Measurement (' num2str(foc,'%.1f') 'mm)'],['Simulation (' num2str(handles.sim_Rcurv(simind),'%i') 'mm, RMSE: ' num2str(err,'%.1f') ')']},'Location','South');
title('Sensitivity (main beam)');

axes(handles.ax_sing_simerr);
hold off;
x= [handles.smod.serr(sensor,:); handles.smod.fielderr(sensor,:)*100];
bar(handles.sim_Rcurv, x');
set(gca,'FontSize',7);
% bar(handles.sim_Rcurv, handles.smod.serr(sensor,:),'EdgeColor','none');
% hold on;
% bar(handles.sim_Rcurv, handles.smod.fielderr(sensor,:),'EdgeColor','none');
xlabel('Curvature (mm)');
ylabel('Deviation from Simulation');


return;
f= figure;
pp = handles.sposa(pos,:);
for j = 1:numel(sens),
    subplot(1,2,1);
    plot(sdist(j)*1e3,sens(j),'.'); hold on;
    xlim([30 50]); ylim([0 1.2]);
    subplot(1,2,2);
    plot(pp(j,1),pp(j,2),'.');hold all;
    xlim([-0.02 0.02]); ylim([-0.02 0.02]);
    pause;
end
close(f);


function handles = plot_sing_sens(handles)
sensor = get(handles.sl_sensor,'Value');
% par.n = str2num(get(handles.e_pixel,'String'));
% par.roi = str2num(get(handles.e_roi,'String'))*1e-3;
par.roi = 40e-3;
par.n = 400;
sens_all = find_sensfield(handles,sensor);

[overlay c2 c3 c6] = find_overlay(handles,sens_all,sensor);

sens_all = round_image(sens_all,14);

% plot
axes(handles.ax_sing_sens); hold off;
image(ind2rgb(round(sens_all*64),jet)); hold all;
h = image(ind2rgb(round(c2*64),gray));
set(h,'AlphaData',c2);
h = image(ind2rgb(round(c3*64),gray));
set(h,'AlphaData',c3);
h = image(ind2rgb(round(c6*64),gray));
set(h,'AlphaData',c6);
axis off;
text(5,15,sprintf('Mismatch: %.2f',overlay*100),'Color','white','FontWeight','bold');

axis on;
set(gca,'FontSize',8);
cent = round(par.n / 2);
tickdst = 0.5e-2;
pertick = par.n / par.roi  * tickdst;
ntick = floor(cent / pertick);
tick = cent-ntick*pertick:pertick:cent+ntick*pertick;
set(gca,'XTick',tick);
set(gca,'XTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
set(gca,'YTick',tick);
set(gca,'YTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
grid on;
title('Sensitivity Field');



function handles = plot_sing_fwhm(handles)
sensor = get(handles.sl_sensor,'Value');
par.roi = 40e-3;
par.n = 400;

fwhm = abs(handles.smod.fwhm(sensor,:));
nax = sqrt(size(handles.sposa,1));
fwhm = reshape(fwhm,nax,nax);

sax = reshape(handles.sposa(:,1),nax,nax);
say = reshape(handles.sposa(:,2),nax,nax);
[X Y] = meshgrid(mean(sax),mean(say'));
dimvec2 = linspace(-par.roi/2,par.roi/2,par.n);
[Xi Yi] = meshgrid(dimvec2,dimvec2);
fwhm_all = interp2(X,Y,fwhm,Xi,Yi,'spline',0);

axes(handles.ax_sing_fwhm);
fwhm_all = round_image(fwhm_all,14,[0.9 3]);
imagesc(fwhm_all);
colormap jet;
caxis([0.9 3]);
axis on;
set(gca,'FontSize',8);
cent = round(par.n / 2);
tickdst = 0.5e-2;
pertick = par.n / par.roi  * tickdst;
ntick = floor(cent / pertick);
tick = cent-ntick*pertick:pertick:cent+ntick*pertick;
set(gca,'XTick',tick);
set(gca,'XTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
set(gca,'YTick',tick);
set(gca,'YTickLabel',num2str((-ntick*tickdst:tickdst:ntick*tickdst)'*1e2,'%.1f'));
grid on;
title('Slice Thickness');


par.cbarfontsize = 7;
par.cbarlabel={'0.9mm','3.0mm'};
plotColorbar(handles.ax_sing_fwhm_cmap,jet,par)


function fwhm_all = round_image(fwhm_all,levels,varargin);
if numel(varargin) >= 1
    x = varargin{1};
    minf = x(1);
    maxf = x(2);
else
    minf = min(fwhm_all(:));
    maxf = max(fwhm_all(:) - minf);
end
fwhm_all = fwhm_all - minf;
fwhm_all(fwhm_all < 0) = 0;
fwhm_all = fwhm_all ./ maxf; 
fwhm_all(fwhm_all > 1) = 1;
fwhm_all = (round(fwhm_all * levels)/levels*maxf) + minf;



function sl_sensor_Callback(hObject, eventdata, handles)
set(handles.sl_sensor,'Value',round(get(handles.sl_sensor,'Value')));
set(handles.e_sensor,'String',num2str(get(handles.sl_sensor,'Value')));
plot_single_sensor(handles,hObject);

function e_sensor_Callback(hObject, eventdata, handles)
set(handles.sl_sensor,'Value',str2num(get(hObject,'String')));
plot_single_sensor(handles,hObject);

function sl_sensor_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function e_sensor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function b_saveindiv_Callback(hObject, eventdata, handles)
dumpAxes(handles.p_ind,handles,sprintf('Sensor_%03i.png',get(handles.sl_sensor,'Value')));
