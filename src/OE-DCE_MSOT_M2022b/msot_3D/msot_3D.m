function varargout = msot_3D(varargin)
% MSOT_3D MATLAB code for msot_3D.fig
%      MSOT_3D, by itself, creates a new MSOT_3D or raises the existing
%      singleton*.
%
%      H = MSOT_3D returns the handle to a new MSOT_3D or the handle to
%      the existing singleton*.
%
%      MSOT_3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSOT_3D.M with the given input arguments.
%
%      MSOT_3D('Property','Value',...) creates a new MSOT_3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msot_3D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msot_3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msot_3D

% Last Modified by GUIDE v2.5 19-Oct-2013 09:38:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msot_3D_OpeningFcn, ...
                   'gui_OutputFcn',  @msot_3D_OutputFcn, ...
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


function msot_3D_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.datainfo = [];
handles.reconWL = 0;
handles.reconSel = 1;
handles.mspSel = 1;
handles.Recon = [];
handles.MSP = [];
handles.grid = false;
handles.mip = 0;
handles.mouse.pressed = false;
handles.bgcmap = gray;
handles.par = struct;
guidata(hObject, handles);

% initialise black axes
axes(handles.ax_3D);
h=image(zeros(100,100,3));
axis off;
set(h,'ButtonDownFcn',@ax_3D_ButtonDownFcn)

axes(handles.ax_bgcbar);
image(zeros(300,100,3));
axis off;

axes(handles.ax_fgcbar);
image(zeros(300,100,3));
axis off;

% parameters
if numel(varargin) >= 1,
   handles.Recon = varargin{1};
   handles.MSP = varargin{2};
   handles.par = varargin{3};
   guidata(hObject, handles);
   
   updateThresholds(handles);
   
   if (isfield(handles.par,'bgthres')),
        set(handles.s_bgthres1,'Value',handles.par.bgthres(1));
        set(handles.s_bgthres2,'Value',handles.par.bgthres(2));
   end
   if (isfield(handles.par,'bgathres')),
        set(handles.s_bgathres1,'Value',round(handles.par.bgathres(1)*64));
        set(handles.s_bgathres2,'Value',round(handles.par.bgathres(2)*64));
   end
   if (isfield(handles.par,'fgthres')),
        set(handles.s_fgthres1,'Value',handles.par.fgthres(1));
        set(handles.s_fgthres2,'Value',handles.par.fgthres(2));
   end
   if (isfield(handles.par,'fgathres')),
        set(handles.s_fgathres1,'Value',round(handles.par.fgathres(1)*64));
        set(handles.s_fgathres2,'Value',round(handles.par.fgathres(2)*64));
   end
   if (isfield(handles.par,'fgmap')),
      cmlist = get(handles.lb_fgcmap,'String');
      cmsel = 0;
      for j = 1:numel(cmlist),
          if strcmp(cmlist{j},handles.par.fgmap),
              cmsel = j;
              break;
          end
      end
      set(handles.lb_fgcmap,'Value',cmsel);
   end
   if (isfield(handles.par,'bgmap')),
      cmlist = get(handles.lb_bgcmap,'String');
      cmsel = 0;
      for j = 1:numel(cmlist),
          if strcmp(cmlist{j},handles.par.bgmap),
              cmsel = j;
              break;
          end
      end
      set(handles.lb_bgcmap,'Value',cmsel);
   end
   if (isfield(handles.par,'zoom')),
       set(handles.s_zoom,'Value',handles.par.zoom);
   end
   if ~isempty(handles.MSP),
       set(handles.cb_fg,'Value',1);
   end
   if isfield(handles.par,'view'),
        set(handles.s_rot1,'Value',handles.par.view(1));
        set(handles.s_rot2,'Value',handles.par.view(2));
        set(handles.s_rot3,'Value',handles.par.view(3));
   end
   
   updateImage(handles);
end


function varargout = msot_3D_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
pos = get(hObject,'Position');
% adapt  anim width
pos_anim = get(handles.p_anim,'Position');
set(handles.p_anim,'Position',[0 0 pos(3) pos_anim(4)]);
% adapt  background width
pos_bg = get(handles.p_bg,'Position');
set(handles.p_bg,'Position',[0 pos_anim(4) pos(3) pos_bg(4)]);
% adapt  foreground width
pos_fg = get(handles.p_fg,'Position');
set(handles.p_fg,'Position',[0 pos_anim(4)+pos_bg(4) pos(3) pos_fg(4)]);
% apat rotation
pos_rot = get(handles.p_rot,'Position');
set(handles.p_rot,'Position',[0 pos_anim(4)+pos_bg(4)+pos_fg(4) pos(3) pos_rot(4)]);
% adapt axes size
pos_ax = get(handles.ax_3D,'Position');
pos_ax(4) = pos(4) - pos_anim(4) - pos_bg(4) - pos_fg(4) - pos_rot(4);
pos_ax(3) = pos(3);
pos_ax(1) = 0;
pos_ax(2) = pos_anim(4)+pos_bg(4)+pos_fg(4)+pos_rot(4);
set(handles.ax_3D,'Position',pos_ax);
updateImage(handles);








% *************************************************************************
% Data Loading
% *************************************************************************


% --------------------------------------------------------------------
function b_opendata_ClickedCallback(hObject, eventdata, handles)
[a b c d] = msot_3D_dataloader(handles.datainfo,...
    handles.reconSel,handles.reconWL,handles.mspSel);
if (isempty(a)) return; end
handles.datainfo = a;
handles.reconSel = b;
handles.reconWL = c;
handles.mspSel = d;

% Load Recon
handles.Recon = loadMSOTRecon(handles.datainfo,handles.reconSel,...
    handles.datainfo.ReconNode(handles.reconSel).ReconStructure(:,:,:,handles.reconWL));

% Load MSP
if (handles.mspSel > 0)
    handles.MSP = loadMSOTMsp(handles.datainfo,handles.mspSel);
    set(handles.lb_fgcomp,'String',handles.datainfo.MSPNode(handles.mspSel).InputSpectra);
    set(handles.lb_fgcomp,'Value',numel(handles.datainfo.MSPNode(handles.mspSel).InputSpectra));
end

guidata(hObject,handles);

updateThresholds(handles);
updateImage(handles);


function updateThresholds(handles)
set(handles.s_bgthres1,'Min',min(handles.Recon(:)));
set(handles.s_bgthres1,'Max',max(handles.Recon(:)));
set(handles.s_bgthres1,'Value',0);
set(handles.s_bgthres2,'Min',min(handles.Recon(:)));
set(handles.s_bgthres2,'Max',max(handles.Recon(:)));
set(handles.s_bgthres2,'Value',max(handles.Recon(:))*0.75);
if (~isempty(handles.MSP))
    cmp = get(handles.lb_fgcomp,'Value');
    msp = handles.MSP(:,:,:,:,:,cmp);
    set(handles.s_fgthres1,'Min',min(msp(:)));
    set(handles.s_fgthres1,'Max',max(msp(:)));
    set(handles.s_fgthres1,'Value',0);
    set(handles.s_fgthres2,'Min',min(msp(:)));
    set(handles.s_fgthres2,'Max',max(msp(:)));
    set(handles.s_fgthres2,'Value',max(msp(:))*0.75);
end

function updateSliders(handles)
set(handles.l_bgthres1,'String',num2str(get(handles.s_bgthres1,'Value'),'%.1e'));
set(handles.l_bgthres2,'String',num2str(get(handles.s_bgthres2,'Value'),'%.1e'));
set(handles.l_bgathres1,'String',num2str(round(get(handles.s_bgathres1,'Value')),'%i'));
set(handles.l_bgathres2,'String',num2str(round(get(handles.s_bgathres2,'Value')),'%i'));
set(handles.l_bgapeak,'String',num2str(round(get(handles.s_bgapeak,'Value')),'%i'));
set(handles.l_bgweight,'String',num2str(get(handles.s_bgweight,'Value'),'%.2f'));
if (~isempty(handles.MSP))
    set(handles.l_fgthres1,'String',num2str(get(handles.s_fgthres1,'Value'),'%.1e'));
    set(handles.l_fgthres2,'String',num2str(get(handles.s_fgthres2,'Value'),'%.1e'));
    set(handles.l_fgathres1,'String',num2str(round(get(handles.s_fgathres1,'Value')),'%i'));
    set(handles.l_fgathres2,'String',num2str(round(get(handles.s_fgathres2,'Value')),'%i'));
    set(handles.l_fgapeak,'String',num2str(round(get(handles.s_fgapeak,'Value')),'%i'));
end
set(handles.l_rot1,'String',num2str(round(get(handles.s_rot1,'Value')),'%i'));
set(handles.l_rot2,'String',num2str(round(get(handles.s_rot2,'Value')),'%i'));
set(handles.l_rot3,'String',num2str(round(get(handles.s_rot3,'Value')),'%i'));
set(handles.l_zoom,'String',num2str(get(handles.s_zoom,'Value'),'%.1f'));
if ~isempty(handles.datainfo),
    set(handles.s_timepoint,'Max',numel(handles.datainfo.Timestamps));
    if numel(handles.datainfo.Timestamps) > 1
        set(handles.s_timepoint,'SliderStep',[1/(numel(handles.datainfo.Timestamps)-1) 1/(numel(handles.datainfo.Timestamps)-1)*5]);
    end
    set(handles.l_timepoint,'String',num2str(handles.datainfo.Timestamps(get(handles.s_timepoint,'Value'))/60,'%.1f min'));
else
    set(handles.s_timepoint,'Enable','off');
end








% *************************************************************************
% Updating the Image
% *************************************************************************

function updateImage(handles,varargin)
updateSliders(handles);
% handles.bgcmap = ones(64,3);
% handles.bgcmap = gray;
if (isempty(handles.Recon))
    return;
end

if (~get(handles.cb_fg,'Value') && ~get(handles.cb_bg,'Value'))
    return;
end

mode = 0;
if (numel(varargin) >= 1)
    mode = varargin{1};
end

[y] = renderImage(handles,mode);
axes(handles.ax_3D);
h=image(y);
axis image;
axis off;
set(h,'ButtonDownFcn',@ax_3D_ButtonDownFcn)



% *************************************************************************
% Rendering the image
% *************************************************************************


function y = renderImage(handles,mode,varargin)

% for animation mode, angle step
addrot = [];
if(numel(varargin) >= 1)
    addrot = varargin{1};
end

fgthres(1) = get(handles.s_fgthres1,'Value');
fgthres(2) = get(handles.s_fgthres2,'Value');

bgthres(1) = get(handles.s_bgthres1,'Value');
bgthres(2) = get(handles.s_bgthres2,'Value');

rotvec(1) = get(handles.s_rot1,'Value');
rotvec(2) = get(handles.s_rot2,'Value');
rotvec(3) = get(handles.s_rot3,'Value');
zoom = get(handles.s_zoom,'Value');
timepoint = get(handles.s_timepoint,'Value');

if ~isempty(handles.datainfo),
    xpx = handles.datainfo.ReconNode(handles.reconSel).ROI / handles.datainfo.ReconNode(handles.reconSel).Resolution*1e3;
    zpx = mean(diff(handles.datainfo.ZPositions));
else
    if isfield(handles.par,'aspect'),
        xpx = handles.par.aspect(1);
        zpx = handles.par.aspect(3);
    else
        xpx = 1;
        zpx = 1;
    end
end

if ~isempty(addrot)
    rotvec = rotvec + addrot;
    rotvec(rotvec > 180) = rotvec(rotvec > 180)-360;
end

clear opt;
opt.ViewerVector = [0 0 1];
opt.Mview=makeViewMatrix(rotvec,[1 1 1/(zpx/xpx)]*zoom,[0 0 0]);
amapsteps = 64;


if (get(handles.cb_bg,'Value'))
    cmlist = get(handles.lb_bgcmap,'String');
%     bgcmap = eval([cmlist{get(handles.lb_bgcmap,'Value')} '(64)']);
    bgcmap = getcmap(cmlist{get(handles.lb_bgcmap,'Value')});
    if ndims(handles.Recon == 3),
        X = handles.Recon;
    else
        X = squeeze(handles.Recon(:,:,timepoint,:,1,1));
    end
    X(X < bgthres(1)) = bgthres(1);
    X(X > bgthres(2)) = bgthres(2);
    X = (X - bgthres(1)) ./ (bgthres(2)- bgthres(1));
    X(X > 0.9999) = 0.9999; % bugfix to prevent hiding of overlay pixels
else
    bgcmap = white;
    X = zeros(size(squeeze(handles.Recon(:,:,timepoint,:,1))));
    bgweight = 1;
end
opt.ColorTable = [bgcmap];

bgalpha = get(handles.s_bgapeak,'Value')/100;
bgalphastart = round(get(handles.s_bgathres1,'Value'));
bgalphaend = round(get(handles.s_bgathres2,'Value'));
bgweight = get(handles.s_bgweight,'Value');
atbg = [zeros(1,bgalphastart), (logspace(0,1,bgalphaend-bgalphastart)-1)/9*bgalpha , ones(1,amapsteps-bgalphaend)*bgalpha];
opt.AlphaTable = [atbg ] ;

if (get(handles.cb_fg,'Value'))
    
    cmp = get(handles.lb_fgcomp,'Value');
    if ndims(handles.MSP) == 3,
        Y = handles.MSP;
    else
        Y = squeeze(handles.MSP(:,:,timepoint,:,1,cmp));
    end
    Y(Y < fgthres(1)) = fgthres(1);
    Y(Y > fgthres(2)) = fgthres(2);
    Y = (Y - fgthres(1)) ./ (fgthres(2)- fgthres(1));
    
    cmlist = get(handles.lb_fgcmap,'String');
%     fgcmap = eval([cmlist{get(handles.lb_fgcmap,'Value')} '(64)']);
    fgcmap = getcmap(cmlist{get(handles.lb_fgcmap,'Value')});
    opt.ColorTable = [bgcmap; fgcmap];    
    
    fgalpha = get(handles.s_fgapeak,'Value')/100;
    fgalphastart = round(get(handles.s_fgathres1,'Value'));
    fgalphaend = round(get(handles.s_fgathres2,'Value'));
    atfg = [zeros(1,fgalphastart), (logspace(0,1,fgalphaend-fgalphastart)-1)/9*fgalpha , ones(1,amapsteps-fgalphaend)*fgalpha];
    opt.AlphaTable = [atbg atfg] ;
    
    % coregister stacks
    if (size(X,3) ~= size(Y,3))
       zpX = handles.datainfo.ZPositions;
       zpY = handles.datainfo.MSPNode(handles.mspSel).ZPositions;
       
       orig = X;
       X = zeros(size(Y));
       for j = 1:size(Y,3)
           ind = find(zpY(j) == zpX);
           if ~isempty(ind)
               X(:,:,j) = orig(:,:,ind(1));
           end
       end
       
    end
end


if handles.mip
    opt.RenderType = 'mip';
else
%     if (get(handles.cb_fg,'Value'))
        opt.RenderType = 'color';
%     else
%         opt.RenderType = 'bw';
%     end
end

% during drag & drop always use fast rendering
if mode == 1, opt.RenderType = 'mip'; end

set(handles.figure1,'Units','pixel');
pos = get(handles.figure1,'Position');
set(handles.figure1,'Units','characters');
n = min(pos(3:4));
opt.ImageSize = [n n];
% opt.ImageSize = [400 400];
opt.WarpInterp = 'nearest';
opt.ShearInterp = 'nearest';

G = zeros(size(X));
if handles.grid
    G = addgrid(G);
    opt.ColorTable = [opt.ColorTable; gray];
    opt.AlphaTable = [opt.AlphaTable linspace(0,1,amapsteps)];
end

if (get(handles.cb_fg,'Value'))
    [Z I] = max([reshape(double(X)*bgweight,prod(size(X)),1) reshape(double(Y),prod(size(Y)),1) reshape(double(G*0.1),prod(size(G)),1)],[],2);
    Z(I == 3) = Z(I == 3)*10 + 2;
    Z(I == 2) = Z(I == 2) + 1;
    Z(I == 1) = Z(I == 1)/bgweight;
    Z = reshape(Z,size(X));
    Z(3,3,1) = 0;
else
    % only bg
    Z = X;
end
% y = render(Z,opt);
warning('off', 'MATLAB:maxNumCompThreads:Deprecated')
% tic
y = render(Z,opt);
% toc
y(y > 1) = 1; y(y<0) = 0;
% size(y)
% max(y(:))


axes(handles.ax_3D);
if (numel(size(y)) == 2)
    y = ind2rgb(round(y*255),gray(255));
end

% plot Colorbars only if real rendermode (mode = 0)
if ~mode
    opt.cbmode = 2;
    opt.cbrotate = 1;
    plotColorbar(handles.ax_bgcbar,[bgcmap atbg'./bgalpha],opt);
    if (get(handles.cb_fg,'Value'))
        plotColorbar(handles.ax_fgcbar,[fgcmap atfg'./fgalpha],opt);
    end
end


% Add outside
function X = addgrid(X)
gv = 0.999999;
gw = 2;
X(1:gw,1:gw,:) = gv;
X(1:gw,end-gw+1:end,:) = gv;
X(end-gw+1:end,end-gw+1:end,:) = gv;
X(end-gw+1:end,1:gw,:) = gv;
X(1:gw,:,1) = gv;
X(end-gw+1:end,:,1) = gv;
X(:,1:gw,1) = gv;
X(:,end-gw+1:end,1) = gv;
X(1:gw,:,end-gw+1:end) = gv;
X(end-gw+1:end,:,end-gw+1:end) = gv;
X(:,1:gw,end-gw+1:end) = gv;
X(:,end-gw+1:end,end-gw+1:end) = gv;














% *************************************************************************
% Background Manipulation
% *************************************************************************

function s_bgthres1_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_bgthres2_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_bgathres1_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_bgathres2_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_bgapeak_Callback(hObject, eventdata, handles)
updateImage(handles);

function cb_bg_Callback(hObject, eventdata, handles)
updateImage(handles);

function lb_bgcmap_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_bgweight_Callback(hObject, eventdata, handles)
updateImage(handles);

function lb_bgcmap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_bgweight_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function s_bgathres1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_bgathres2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_bgapeak_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_bgthres1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_bgthres2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end








% *************************************************************************
% Foreground Manipulation
% *************************************************************************

function s_fgthres1_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_fgthres2_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_fgathres1_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_fgathres2_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_fgapeak_Callback(hObject, eventdata, handles)
updateImage(handles);

function cb_fg_Callback(hObject, eventdata, handles)
updateImage(handles);

function lb_fgcomp_Callback(hObject, eventdata, handles)
updateThresholds(handles);
updateImage(handles);

function lb_fgcmap_Callback(hObject, eventdata, handles)
updateImage(handles);


function s_fgthres1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_fgthres2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_fgathres1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_fgathres2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_fgapeak_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function lb_fgcomp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb_fgcmap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end














% *************************************************************************
% View Manipulation
% *************************************************************************

function s_rot1_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_rot2_Callback(hObject, eventdata, handles)
updateImage(handles);

function s_rot3_Callback(hObject, eventdata, handles)
updateImage(handles);


function s_rot1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_rot2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function s_rot3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function l_rot1_Callback(hObject, eventdata, handles)
set(handles.s_rot1,'Value',str2num(get(hObject,'String')));
updateImage(handles);

function l_rot2_Callback(hObject, eventdata, handles)
set(handles.s_rot2,'Value',str2num(get(hObject,'String')));
updateImage(handles);

function l_rot3_Callback(hObject, eventdata, handles)
set(handles.s_rot3,'Value',str2num(get(hObject,'String')));
updateImage(handles);


function s_zoom_Callback(hObject, eventdata, handles)
updateImage(handles);

function l_zoom_Callback(hObject, eventdata, handles)
updateImage(handles);


function s_zoom_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function l_zoom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s_timepoint_Callback(hObject, eventdata, handles)
set(hObject,'Value',round(get(hObject,'Value')));
updateImage(handles);

function s_timepoint_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% -------------------------------------------------------------------------
% Drag & Drop
% -------------------------------------------------------------------------

% Start dragging
function ax_3D_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ax_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(get(hObject,'Parent'));
handles.mouse.pressed = true;
handles.mouse.oldpos = get(handles.figure1,'CurrentPoint');
% get(get(hObject,'Parent'),'CurrentPoint')
guidata(get(hObject,'Parent'),handles);


% Dragging finished
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
updateAngles(hObject,handles,0);

% During dragging
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
updateAngles(hObject,handles,1);


% Calculate new angles after drag and drop
function updateAngles(hObject,handles,cc)
if (handles.mouse.pressed)
    newpos = get(hObject,'CurrentPoint');
    pdiff = newpos - handles.mouse.oldpos;
    rev = 1;
    if (abs(get(handles.s_rot1,'Value')) > 90) rev = -1; end
    newzrot = get(handles.s_rot1,'Value') + atan(pdiff(1)/50)*360/pi*rev;
    if(newzrot > 180) newzrot = newzrot - 360; end
    if(newzrot < -180) newzrot = newzrot + 360; end
    set(handles.s_rot1,'Value',newzrot);

    rev = 1;
    if (abs(get(handles.s_rot2,'Value')) > 90) rev = -1; end
    newyrot = get(handles.s_rot2,'Value') + atan(pdiff(2)/50)*360/pi*rev;
    if(newyrot > 180)   newyrot = newyrot - 360;    end
    if(newyrot < -180)  newyrot = newyrot + 360;    end
    set(handles.s_rot2,'Value',newyrot);

    if cc == 0
        handles.mouse.pressed = false;
        guidata(hObject,handles);
    end
  
    updateImage(handles,cc);
end


% -------------------------------------------------------------------------
% Toolbar Buttons
% -------------------------------------------------------------------------
function t_mip_ClickedCallback(hObject, eventdata, handles)
if (strcmp(get(hObject,'State'),'on'))
    handles.mip = true;
else
    handles.mip = false;
end
guidata(hObject,handles);
updateImage(handles);

function t_grid_ClickedCallback(hObject, eventdata, handles)
if (strcmp(get(hObject,'State'),'on'))
    handles.grid = true;
else
    handles.grid = false;
end
guidata(hObject,handles);
updateImage(handles);



% -------------------------------------------------------------------------
% Animation
% -------------------------------------------------------------------------


function e_nframes_Callback(hObject, eventdata, handles)
function e_nframes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_step3_Callback(hObject, eventdata, handles)
function e_step3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_step2_Callback(hObject, eventdata, handles)
function e_step2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function e_step1_Callback(hObject, eventdata, handles)
function e_step1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function b_prevanim_Callback(hObject, eventdata, handles)

nframes = str2num(get(handles.e_nframes,'String'))
step(1) = str2num(get(handles.e_step1,'String'));
step(2) = str2num(get(handles.e_step2,'String'));
step(3) = str2num(get(handles.e_step3,'String'))

wbar = waitbar(0,'Rendering Images');
for j = 1:nframes
    waitbar(j/nframes,wbar);
    img = renderImage(handles,2,(j-1)*step);
    y(:,:,:,j) = img;
end
close(wbar);


f = figure;
set(f,'Units','pixel');
pos = get(f,'Position');
set(f,'Position',[200 200 size(y,1) size(y,2)]);
whitebg(f);
ax = axes('Units','pixel','Position',[1 1 size(y,1) size(y,2)]);
c = 0;
while(ishandle(f))
    c = c + 1;
    if c > nframes, c = 1; end
    axes(ax),image(y(:,:,:,c));
    pause(0.15);
end


function b_saveanim_Callback(hObject, eventdata, handles)
nframes = str2num(get(handles.e_nframes,'String'))
step(1) = str2num(get(handles.e_step1,'String'));
step(2) = str2num(get(handles.e_step2,'String'));
step(3) = str2num(get(handles.e_step3,'String'))


[fname p] = uiputfile('*.tiff')
if isempty(fname)
    return;
end


wbar = waitbar(0,'Rendering Images');
for j = 1:nframes
    waitbar(j/nframes,wbar);
    img = renderImage(handles,2,(j-1)*step);
    y(:,:,:,j) = img;
end
close(wbar);

if(exist([p fname],'file'))
    delete([p fname]);
end

wbar = waitbar(0,'Saving Images');
for j = 1:nframes
    waitbar(j/nframes,wbar);
    imwrite(y(:,:,:,j),[p fname],'WriteMode','append');
end
close(wbar);

