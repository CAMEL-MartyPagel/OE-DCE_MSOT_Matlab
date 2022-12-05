function varargout = daq_chkframes(varargin)
% DAQ_CHKFRAMES MATLAB code for daq_chkframes.fig
%      DAQ_CHKFRAMES, by itself, creates a new DAQ_CHKFRAMES or raises the existing
%      singleton*.
%
%      H = DAQ_CHKFRAMES returns the handle to a new DAQ_CHKFRAMES or the handle to
%      the existing singleton*.
%
%      DAQ_CHKFRAMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAQ_CHKFRAMES.M with the given input arguments.
%
%      DAQ_CHKFRAMES('Property','Value',...) creates a new DAQ_CHKFRAMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before daq_chkframes_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to daq_chkframes_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help daq_chkframes

% Last Modified by GUIDE v2.5 16-Oct-2014 10:04:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @daq_chkframes_OpeningFcn, ...
                   'gui_OutputFcn',  @daq_chkframes_OutputFcn, ...
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


% --- Executes just before daq_chkframes is made visible.
function daq_chkframes_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to daq_chkframes (see VARARGIN)

% Choose default command line output for daq_chkframes
handles.output = hObject;
handles.indir = '';
handles.logfid = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes daq_chkframes wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = daq_chkframes_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function b_opendir_Callback(hObject, eventdata, handles)
newpath = uigetdir(handles.indir);
if (~exist(newpath,'dir'))
    errordlg('cannot find folder');
    return;
end

handles.indir = newpath;
guidata(hObject,handles);

b_run_Callback(hObject, eventdata, handles);
return;


% --- Executes on button press in b_run.
function b_run_Callback(hObject, eventdata, handles)
if isempty(handles.indir), return; end

set(handles.l_foldername,'String',handles.indir);
set(handles.lb_output,'Value',1);
set(handles.lb_output,'String','');
set(handles.lb_dirlist,'Value',1);
set(handles.lb_dirlist,'String','');
set(handles.lb_direlist,'Value',1);
set(handles.lb_direlist,'String','');

handles.logfid = fopen([handles.indir '\daq_check.log'],'a');
fprintf(handles.logfid,'\n\n\n\n\nRunning daq_check...\n');

guidata(hObject,handles);

% run analysis
analyse_daq_frames(handles.indir,handles);

% close logfile
fclose(handles.logfid);
handles.logfid = 0;
guidata(hObject,handles);

    

function lb_output_Callback(hObject, eventdata, handles)


function lb_output_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Analyse DAQ Frames
function analyse_daq_frames(sdir,handles)

sfiles = dir([sdir '\*.bin']);
fnames = {sfiles.name};
[tmp ind] = sort(fnames);
sfiles = sfiles(ind);

dt = 1/str2num(get(handles.e_rep,'String'))*1e3;
set(handles.l_size,'String',sprintf('%i (%.1f min)',numel(sfiles),numel(sfiles)*dt/60000));
% datearr = cell2mat({sfiles.datenum});
% datearr = datearr-min(datearr(:))
%
fastmode = get(handles.cb_fast,'Value');
if fastmode
    logline('Operating in fast mode...\n',handles);
end


wbar = waitbar(0,'Reading Frame Files','Name','Analyzing','CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
lastseq = 0;
lastts = 0;
seqvec = zeros(1,numel(sfiles));
seqvec(:) = nan;
flist = cell(numel(sfiles),1);
felist = cell(0,0);
for sfid = 1:numel(sfiles)
    
    if getappdata(wbar,'canceling')
        logline(sprintf('Analsis cancelled at number %i\n',sfid),handles);
        errordlg(sprintf('Analysis cancelled at number %i',sfid));
        break;
    end
    
    nchannels = sfiles(sfid).bytes/2032/2;
    
    FID = fopen([sdir '\' sfiles(sfid).name],'r');
    lbyte1 = 0;
    lbyte2 = 0;
    lbyte3 = 0;
    err = 0;
    
    if fastmode
        % read beginning of frame
        lbyte1 = fread(FID,1,'uint16');
        lbyte2 = fread(FID,1,'uint16');
        lbyte3 = fread(FID,1,'uint16');

        % move to end of frame and read last
        fseek(FID,2032*2*(nchannels-2)+2029*2,0);
        byte1 = fread(FID,1,'uint16');
        byte2 = fread(FID,1,'uint16');
        byte3 = fread(FID,1,'uint16');

        if byte1 ~= lbyte1+nchannels-1
            logline(sprintf('Channel number mismatch: first %i, last %i in seq %i (%s)\n',lbyte1,byte1,byte2,sfiles(sfid).name),handles);
            err = 1;
        end
        if byte2 ~= lbyte2
            logline(sprintf('Sequence number mismatch within frame %s (first: %i last: %i)\n',sfiles(sfid).name,lbyte2,byte2),handles);
            err = 1;
        end
        if get(handles.cb_timestamp,'Value') && lbyte3 && byte3 ~= lbyte3
            logline(sprintf('Timestamp mismatch within frame %s (@%i after @%i)\n',sfiles(sfid).name,byte3,lbyte3),handles);
            err = 1;
        end
    else
    
        for ch = 1:nchannels
            byte1 = fread(FID,1,'uint16');
            byte2 = fread(FID,1,'uint16');
            byte3 = fread(FID,1,'uint16');
            
            if byte1 ~= lbyte1+1
                logline(sprintf('Channel number sorting mismatch: %i after %i in seq %i (%s)\n',byte1,lbyte1,byte2,sfiles(sfid).name),handles);
                err = 1;
            end
            if lbyte2 && byte2 ~= lbyte2
                logline(sprintf('Sequence number mismatch within frame %s at channel %i (%i after %i)\n',sfiles(sfid).name,byte1,byte2,lbyte2),handles);
                err = 1;
            end
            if get(handles.cb_timestamp,'Value') && lbyte3 && byte3 ~= lbyte3
                logline(sprintf('Timestamp mismatch within frame %s at channel %i (@%i after @%i)\n',sfiles(sfid).name,byte1,byte3,lbyte3),handles);
                err = 1;
            end

            % prepare for next
            fseek(FID,2029*2,0);
            lbyte1 = byte1;
            lbyte2 = byte2;
            lbyte3 = byte3;
        end
    end % fast mode
    fclose(FID);
    
    if lastseq && lbyte2 ~= lastseq + 1
      logline(sprintf('Sequence mismatch (%i after %i) in %s\n',lbyte2,lastseq,sfiles(sfid).name),handles);  
      err = 1;
    end
    seqvec(sfid) = lbyte2;
    lastseq = lbyte2;
    
    if get(handles.cb_timestamp,'Value'); 
        if lastts && lbyte3 ~= lastts + dt,
          logline(sprintf('Timestamp mismatch (%i after %i) in %s\n',lbyte3,lastts,sfiles(sfid).name),handles);  
          err = 1;
        end
    else
        lbyte3 = 0; 
    end;
    seqts(sfid) = lbyte3;
    lastts = lbyte3;    
    
    flist{sfid} = addframe(sfid,lastseq,lbyte3,sfiles(sfid).name,err);
    if err, felist{numel(felist)+1} = flist{sfid}; end;

    % progress bar
    wbar = waitbar(sfid/numel(sfiles),wbar,sprintf( 'Analyzing %i of %i', sfid, numel(sfiles) )) ;

end
close(wbar);delete(wbar);
set(handles.lb_dirlist,'String',char(flist{:}));
set(handles.lb_direlist,'String',char(felist{:}));

 t = 0:0.1:0.1*(numel(seqvec)-2);
 sdiff = diff(seqvec);
 tdiff = diff(seqts);
 tdiff(tdiff == -(2^16-dt)) = dt;  % timestamp reset
 tdiff = uint16(tdiff);
figure('Position',[250 100 1550 600]);
[ax h1 h2] = plotyy(t/60,sdiff,t/60,tdiff);
set(h1,'Marker','.','MarkerSize',15,'LineStyle','none');
set(h2,'Marker','x','MarkerSize',10,'LineStyle','none');
xlabel('time (min)');
ylabel('sequence difference');
line(xlim,[0 0]);
ylim(ax(1),[min(sdiff)-1 max(sdiff)+1]);
ylim(ax(2),[min(tdiff)-1 max(tdiff)+1]);
title(sdir,'Interpreter','none');
% logline(num2str(seqvec),handles);



function logline(str,handles)
s = cellstr(get(handles.lb_output,'String'));
s = [s ; cellstr(str)];
set(handles.lb_output,'String',s);
if (handles.logfid) 
    fprintf(handles.logfid,str);
end

function s = addframe(nb,seq,ts,str,err)
s = [num2str(nb,'%06i') ' - ' num2str(seq,'%04i')];
if ts > 0, s = [s '@' num2str(ts,'%05i')]; end;
s = [s ': ' str];
if err
    s = ['* ' s];
end


function cb_fast_Callback(hObject, eventdata, handles)
function lb_dirlist_Callback(hObject, eventdata, handles)
function lb_dirlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function lb_direlist_Callback(hObject, eventdata, handles)
function lb_direlist_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function cb_timestamp_Callback(hObject, eventdata, handles)


function cb_filter_Callback(hObject, eventdata, handles)
if (get(hObject,'Value')),
    set(handles.lb_direlist,'Visible','on');
    set(handles.lb_dirlist,'Visible','off');
else
    set(handles.lb_dirlist,'Visible','on');
    set(handles.lb_direlist,'Visible','off');
end

function b_showframe_Callback(hObject, eventdata, handles)
if get(handles.cb_filter,'Value');
    str = get(handles.lb_direlist,'String');
    select = get(handles.lb_direlist,'Value');
else
    str = get(handles.lb_dirlist,'String');
    select = get(handles.lb_dirlist,'Value');
end
str = str(select,:);
sep = strfind(str,': ');
fname = str(sep(1)+2:end);

d = dir([handles.indir '\' fname]);
fsize = d.bytes;
nchannels = fsize/2032/2;
FID = fopen([handles.indir '\' fname],'r');

if get(handles.cb_timestamp,'Value'),
    fsize = 2029;
else
    fsize = 2030;
end

dat = zeros(fsize,nchannels);
for ch = 1:nchannels
    byte1(ch) = fread(FID,1,'uint16');
    byte2(ch) = fread(FID,1,'uint16');
    if get(handles.cb_timestamp,'Value'),
        byte3(ch) = fread(FID,1,'uint16');
    else
        byte3(ch) = 0;
    end

    % prepare for next
    dat(:,ch) = fread(FID,fsize,'uint16');
end
fclose(FID);

[b_HPF,a_HPF] = cheby1( 1, 1, 2 * 500e3/4e7 , 'high' ) ; 
datf = FilterM( b_HPF, a_HPF, dat );

f=figure;
subplot(1,2,1);
imagesc(datf);
% colorbar;
title(str,'interpreter','none');

subplot(2,2,2);
plotyy(1:numel(byte1),byte1,1:numel(byte1)-1,diff(byte1));
title('channel numbers (byte 1)');
% xlim([0 256]);

subplot(2,2,4);
[ax h1 h2] = plotyy(1:numel(byte2),byte2,1:numel(byte3),byte3);
ylim(ax(1),[min(byte2)-1 max(byte2)+1]);
ylim(ax(2),[min(byte3)-1 max(byte3)+2]);
n=unique(uint16(get(ax(2),'ytick')));
set(ax(2),'ytick',n);
set(ax(2),'yticklabel',num2str(n','%i'));

% xlim([0 numel(byte1)]);
title('sequence number');
legend({'sequence (byte2)','timestamp (byte3)'});



function e_rep_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function e_rep_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
