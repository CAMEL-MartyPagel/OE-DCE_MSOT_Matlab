function [R, wls, zpos, ts, sigMatF, sigMat, par, spar] = reconMSOT2(fn,varargin)
% [R, wls, zpos, ts, sigFilt, sig, parameters] = reconMSOT2(INPUT[,parameters,verbose]);
% 
% Reconstruction function to reconstruct a selected dataset using the
% supplied parameters using the object oriented OpenCL implementation
% NOTE: Only backprojection supported so far
%
% Input Parameters:
% -----------------
% - INPUT: Either 
%      * filename (char or cell) of .msot file 
%      * datainfo object
%      * signal matrix (requires mandatory parameters)
% - Parameters: (struct)
%   image_select -> Reconstruction algorithm (*direct, model_lin, wavelet)
%   n            -> Resolution (default: 250)
%   proj         -> Number of Projections (default: 512)
%   c            -> Speed of Sound (either absolute or offset, default: 0)
%   roi          -> Region of Interest x/y (in m, default: 25)
%   roiZ         -> Region of Interest z (in m, default: 25)
%   center       -> Center of reconstructed volume (in m, default: [0 0 0]
%                     REMARK: Currently only works for z
%   limits       -> Direct input for CL library - overrides ROI Settings
%   pos_sensor   -> Sensor positions as 3D vector (in m)
%   rotate_sensor-> Rotate Sensor in degrees (only for parametric desc)
%   limit_sensors-> Specify sensors to use (to reduce or reorder sensors)
%   fs           -> Sampling Frequency (for filters)
%   filter_f     -> Bandpass frequencies (default: 50kHz - 6.5Mhz)
%   filter_rise  -> Width of rising slope (percent of LCF, default: 0.2)
%   filter_fall  -> Width of falling slope (percent of UCF, default: 0.2)
%   freqVec      -> Frequency Vector as direct input for CL - overrides
%                     Filters and requires ampVec
%   ampVec       -> Amplitude Vector - see above
%   impresp      -> Specify Impulse Response 
%                     0: none, 
%                     []: from scan, 
%                     char: filename of file, 
%                     double array: time domain resp,
%                     3 element vector: Center Frequency, BW, Phase Shift
%   timeres      -> Timeresolution 
%   selMat       -> Selection Matrix based on datainfo.ScanStructure
%   CLDir        -> Directory of itheract.cl (alternatively, specify
%                     environment variable ITHERFILES1)
%   CLVendor     -> OpenCL device vendor: 'intel', 'amd' or 'nvidia'
%   CLType       -> OpenCL device type: 'cpu' or 'gpu'
%   progress     -> Show progress window (default: true)
%   verbose      -> Print verbose logging (default: false)
%   save         -> Save reconstructed images as .mat file (default: false)
%
% Return Values:
% --------------
% @return R          An array containing the reconstructed images in 8D
%                    1) x
%                    2) y
%                    3) z (only if 3D detector, currently unused)
%                    4) RUN (repetition of the whole acquisition step, time
%                       dimension (timepoints as separate vector)
%                    5) z-position (if 2D system with translation stage)
%                    6) REPETITION (currently unused)
%                    7) Wavelength (see wls vector)
%                    8) Individual Images (only if not averaged)
% @return wls        Wavelength Vector
% @return zpos       Array of z-stage positions
% @return ts         Multidimensional array of timestamps in s, same 
%                    dimension as R
% @return sigFilt    Filtered Signals in single precision
% @return sig        Raw Signals, scaled to uint16
% @return parameters Complete Parameter Struct

%% Input Parameters
% Filename, Signal or Datainfo
sigMat = [];
spar = struct;
datainfo = [];

if ( isstruct(fn) || isa(fn,'msotData') )
    datainfo = fn;
    fn = datainfo.XMLFileName;
elseif (isnumeric(fn) || isinteger(fn)),
    sigMat = fn;
else
    if (~exist(fn,'file'))
      error(['MSOT File not found: ' fn]);
    else 
        datainfo = loadMSOT(fn);
    end
end

if numel(varargin) >= 1 && isstruct(varargin{1})
    cpar = varargin{1};
else
    cpar = [];
end

% defaults
par = struct();
par.image_select = 'direct';
par.n = 250;
par.nZ = 1;
par.proj = 512;
par.c = 1510;
par.c_round = false;
par.ua = 0;
par.roi = 25e-3;
par.roiZ = 25e-3;
par.center = [0 0 0];    % center of reconstruction (currently only works for z)
par.limits = [];    % limits vector (overrides the roi Settings)
par.impresp = [];   % will be loaded automatically if included w scan, put
                    % 0 if desired to NOT include IR

par.fs = 0;         % sampling frequency 
par.pos_sensor = [];        % sensor positions (3D)
par.rotate_sensor = 0;      % rotate sensor angles (in degree)
par.limit_sensors = [];  % subset of sensors (can also be used to reorder)
par.sensor_sensitivity = [];  % sensitivity vector to be divided from sensors to even out sensivity differences
par.interpolate_sensors = []; % subset of sensors to interpolate from neighboring channels

par.directbp = false;           % use real direct BP (not derivative)
par.phaseflip = true;           % phase flipping (active in vM 3.8 for KLM model EIR)
par.filter_f = [50 6500]*1e3;
par.filter_rise = 0.2;          % percentage of rise of FIR filter (width of slope)
par.filter_fall = 0.2;          % percentage of rise of FIR filter (width of slope)
par.timeres = 3;
par.freqVec = [];               % manual input for filter vector
par.ampVec = [];                % manual input for amplitude vector

par.lsqr_reg = 1;               % choose regularized model-based

par.mirror = false; % mirror image (0: none, true/3: both axis, 1: x, 2: y)

par.CouplantSoS = 0;
par.CouplantRadius = 0;
par.CouplantOrigin = 0;

par.CLDir = '';
par.CLVendor = '';      % can be 'intel', 'amd' or 'nvidia'
par.CLType = '';        % can be 'cpu' or 'gpu'

% slightly change defaults if datainfo is present
if ~isempty(datainfo),
    par.selMat = datainfo.ScanStructure;         % load whole reconstruction
    par.fs = datainfo.HWDesc.SamplingFrequency;

    
    % default speed of sound from temperature
    T = datainfo.AverageTemperature;
    par.c = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;

    % if c offset is given in parameters make sure it gets used
    if ~isempty(cpar) && isfield(cpar,'c') && cpar.c < 1000,
        par.c = par.c + cpar.c;
        cpar.c = par.c;
    end
    
    % get sensor positions
    if isempty(par.pos_sensor),
        if isfield(datainfo.HWDesc,'ProjectionData') && ~isempty(datainfo.HWDesc.ProjectionData),
            par.pos_sensor = datainfo.HWDesc.ProjectionData;
            if mean(datainfo.HWDesc.ProjectionData(:,2)) == 0, % new geometry for 4.0
                par.pos_sensor(:,1) = -datainfo.HWDesc.ProjectionData(:,1);
                par.pos_sensor(:,2) = -datainfo.HWDesc.ProjectionData(:,3);
                par.pos_sensor(:,3) = -datainfo.HWDesc.ProjectionData(:,2);
            end
        else
            par.pos_sensor = linspace(datainfo.HWDesc.StartAngle,datainfo.HWDesc.EndAngle,datainfo.HWDesc.NumDetectors)+par.rotate_sensor/360*pi*2;
            par.pos_sensor = [ cos( par.pos_sensor' ) sin( par.pos_sensor' ) zeros( size( par.pos_sensor' ) ) ];
            par.pos_sensor = par.pos_sensor .* ( (ones(size(par.pos_sensor,1),1) .* datainfo.HWDesc.Radius) * ones(1,3) );
        end
        par.center(1) = datainfo.HWDesc.Radius;
    end
    
    % if 3D dataset, automatically do 3D recon
    if datainfo.is3D, par.nZ = par.n; end;
end

% specials
par.save = false;       
par.progress = true;
par.verbose = false;


% Copy parameters from input struct
if ~isempty(cpar),
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

% if signal Matrix is supplied
if isfield(par,'sigMat'),
    sigMat = par.sigMat;
end

% validate
if par.fs == 0, error('Please set sampling frequency in par.fs'); end
if size(par.pos_sensor,2) ~= 3, error('Please supply par.pos_sensor as 3D vector'); end;

%% check GPU environment

% CL File
oldcldir = getenv('ITHERAFILES1'); % backup environment
clfile = which('itheract.cl');     % find on path
% path directly specified
if exist([par.CLDir '\itheractl.cl'],'file') == 2,
    setenv('ITHERAFILES1',par.CLDir);
% using environment variable
elseif ~isempty(getenv('ITHERAFILES1')) && exist([getenv('ITHERAFILES1') '\itheract.cl'],'file') == 2,
    par.CLDir = getenv('ITHERAFILES1');
% maybe its on path?
elseif ~isempty(clfile),
    [clpath,~] = fileparts(clfile);
    par.CLDir = clpath;
% setting environment variables at runtime wont help the itheractcl.dll
% until passing an explicit path to the CL file is possible, work around
% this by copying the cl file to a specific location and issue a warning
%     setenv('ITHERAFILES1',clpath);
%     try
%         copyfile(clfile,getenv('TMP'));
%     catch
%     end
%     fprintf('Warning: CL Sources might not be found. Point ITHERAFILES1 to %s\n',getenv('TMP'));
%     par.CLDir = getenv('TMP');
% otherwise fail
else
    error('itheract.cl is not found in your environment\nPlease use "setenv ITHERAFILES1 YOURPATH" or supply using par.CLDir');
end

if par.verbose, fprintf('Using CL Sources in %s (ITHERAFILES1: %s)\n',par.CLDir,getenv('ITHERAFILES1')); end    


% Find MEX Files
if exist('startupcl','file') ~= 3,
   % check in ITHERAFILES1
  if exist([getenv('ITHERAFILES1') '\startupcl.mexw64'],'file') == 3,
      path(getenv('ITHERAFILES1'),path);
      if par.verbose, fprintf('Automatically using mex-Files in %s\n',oldcldir); end    
  end
  % check again
  if exist('startupcl','file') == 3,
  else
      error('startupcl.mexw64 is not on your path\n');
  end
end
if exist('clBackprojectionSync','file') ~= 3,
  error('clBackprojectionSync.mexw64 is not on your path\n');
end
if exist('cleanupcl','file') ~= 3,
  error('cleanupcl.mexw64 is not on your path\n');
end


% DLL 
if exist('itheractcl.dll','file') ~= 2,
  error('itheractcl.dll is not on your path\n');
end


% MATLAB Files
if exist('clBackprojection','file') ~= 2,
  error('clBackprojection.m is not on your path\n');
end

% CL Options
if isempty(par.CLVendor) && ~isempty(getenv('ITHERA_CL_VENDOR')),
    par.CLVendor = getenv('ITHERA_CL_VENDOR');
end
if isempty(par.CLType) && ~isempty(getenv('ITHERA_CL_TYPE')),
    par.CLType = getenv('ITHERA_CL_TYPE');
end

%% Start CL
if par.verbose, par.CLDebug = '1'; fprintf('Enabling CL Debugging\n'); else par.CLDebug = '0'; end;
if ~isempty(par.CLVendor) && ~isempty(par.CLType),
    if par.verbose, fprintf('Using CL Parameters: vendor %s, type %s\n',par.CLVendor,par.CLType); end
    startupcl('vendor',par.CLVendor,'type',par.CLType,'clpath',par.CLDir,'debuglevel',par.CLDebug);
elseif ~isempty(par.CLVendor)
    if par.verbose, fprintf('Using CL Parameters: vendor %s\n',par.CLVendor); end
    startupcl('vendor',par.CLVendor,'clpath',par.CLDir,'debuglevel',par.CLDebug);
elseif ~isempty(par.CLType),
    if par.verbose, fprintf('Using CL Parameters: type %s\n',par.CLType); end
    startupcl('type',par.CLType,'clpath',par.CLDir,'debuglevel',par.CLDebug);
else
    if par.verbose, fprintf('Using CL Parameters: none\n',par.CLType,'clpath',par.CLDir,'debuglevel',par.CLDebug); end
    startupcl('clpath',par.CLDir,'debuglevel',par.CLDebug);
end


%% load signal Data
if isempty(sigMat),
    [sigMat,~,selMat,spar] = loadMSOTSignals(datainfo,par.selMat,par);
else
    selMat = par.selMat; 
end
% ts = datainfo.RelTime(selMat);

% put selMat to 1D 
sigMat = reshape(sigMat,size(sigMat,1),size(sigMat,2),numel(selMat));
svec = size(selMat);
selMat = reshape(selMat,numel(selMat),1);

% truncate/reorder data and sensor positions 
if ~isempty(par.limit_sensors)
    par.pos_sensor = par.pos_sensor(par.limit_sensors,:);
    sigMat = sigMat(:,par.limit_sensors,:);
    par.proj = numel(par.limit_sensors);
    if par.verbose, fprintf('Signals and Sensor Positions truncated to %i sensors\n',numel(par.limit_sensors)); end;
end

% interpolate faulty sensors
if ~isempty(par.interpolate_sensors)
    fprintf('Interpolating for faulty Sensors...\n');
    loc = ones(size(sigMat,2),1,'logical');
    loc(par.interpolate_sensors) = logical(0);
    sensidx = 1:size(sigMat,2);
    [X_,Y_] = meshgrid(sensidx(loc),1:2030);
    [X,Y] = meshgrid(sensidx,1:2030);
    for jj = 1:size(sigMat,3)
        sigMat(:,:,jj) = interp2(X_,Y_,sigMat(:,loc,jj),X,Y);
    end
end

% amplify sensors
if ~isempty(par.sensor_sensitivity)
    fprintf('Incorporating sensitivity vector for Sensors...\n');
    for jj = 1:size(sigMat,3)
        sigMat(:,:,jj) = sigMat(:,:,jj) ./ ( ones(size(sigMat,1),1) * reshape( par.sensor_sensitivity,1,numel(par.sensor_sensitivity) )  ) ;
    end
end
    

% convert to uint16
if strcmp(par.image_select,'direct') || strcmp(par.image_select,'model_lin') || strcmp(par.image_select,'nonNeg_MB'),
if isinteger(sigMat),
    fac = single(1);
    scalefac = [datainfo.ScanFrames(selMat).CorrectionFactor];
else
    sigMat = sigMat - min(sigMat(:));   % remove negatives
    fac = floor((2.^16)./max(sigMat(:)));
    sigMat = uint16(sigMat * fac);
    scalefac = ones(numel(selMat),1);
end
end

% round SoS
if par.c_round,
    par.c = round(par.c);
end

% supply ua
if numel(par.ua) ~= numel(selMat),
    par.ua = ones(numel(selMat),1)*par.ua;
end

%% load impulse response
if par.impresp == 0, % avoid impresp
    par.impresp = [];
elseif isa(par.impresp,'msotEIR'),
    par.impresp = par.impresp.tdom;
elseif isnumeric(par.impresp) && numel(par.impresp) == 2,
    eir = msotEIR(par.impresp(1),par.impresp(2));
    par.impresp = eir.tdom;
elseif isnumeric(par.impresp) && numel(par.impresp) == 3,
    eir = msotEIR(par.impresp(1),par.impresp(2),par.impresp(3));
    par.impresp = eir.tdom;
elseif isnumeric(par.impresp) && numel(par.impresp) == 4,
    eir = msotEIR(par.impresp(1),par.impresp(2),par.impresp(3),par.impresp(4));
    par.impresp = eir.tdom;
elseif ischar(par.impresp) || isempty(par.impresp),
    if ischar(par.impresp),
        irpath = par.impresp;
    else
        irpath = strrep((datainfo.XMLFileName),'msot','irf');
    end
    
    if exist(irpath,'file'),
        FID = fopen(irpath);
        par.impresp = fread(FID,'double')';
        fclose(FID);
        if par.verbose, fprintf('Electrical Impulse Response loaded (%s)\n',irpath); end
    else 
        par.impresp = [];
    end
end

if isnumeric(par.impresp) && numel(par.impresp) == 2030 && par.impresp(1) >= 1000,
    par.impresp(1) = par.impresp(1) - 1000;
    fprintf('Auto-Detecting Direct Backprojection\n');
    par.directbp = true;
end

if isempty(par.impresp) && par.verbose,
   fprintf('No Electrical Impulse Response used\n');
end    


%% preparation

if strcmp(par.image_select,'direct') || strcmp(par.image_select,'model_lin') || strcmp(par.image_select,'nonNeg_MB'),
    if (isempty(par.limits)),
        par.limits(1) = par.center(1) - par.roi/2;              % x low
        par.limits(2) = par.center(1) + par.roi/2;        % x high
        par.limits(3) = par.center(3) - par.roiZ / 2;    % FoV start
        par.limits(4) = par.center(3) + par.roiZ / 2;    % FOV end
         
    end
    
    % automatically design filter vector based on cut-off frequencies
    if isempty(par.freqVec)
        f1=par.filter_f(1)*1e-6;
        f2=par.filter_f(2)*1e-6;
        if (f2 == 0), f2 = 12; end;

        if par.directbp && strcmp(par.image_select,'direct') && par.phaseflip,
            ampmin = 1.0*1i;
            ampmax = 1.0*1i;
        elseif par.directbp && (strcmp(par.image_select,'model_lin') || strcmp(par.image_select,'nonNeg_MB') || strcmp(par.image_select,'direct')),
            ampmin = 1.0;
            ampmax = 1.0;
        elseif ~par.directbp && strcmp(par.image_select,'direct'),
            ampmin = (1/20*f1)*1i;
            ampmax = (1/20*f2)*1i;
        elseif ~par.directbp && (strcmp(par.image_select,'model_lin') || strcmp(par.image_select,'nonNeg_MB')),
            ampmin = (1/20*f1);
            ampmax = (1/20*f2);
        end
 
        par.freqVec = [0.0, (1-par.filter_rise)*f1, f1,  f2, f2+par.filter_fall*f2, 20]/20;   % frequency
        par.ampVec =  [0.0, 0.0, ampmin, ampmax, 0.0, 0.0];   
    end

    if par.nZ == 1 && par.CouplantSoS == 0,
        type = '2d';
    else
        if par.CouplantRadius ~= 0,  type = '3d_c';     % curved couplant
        elseif par.CouplantSoS ~= 0, type = '3d_f';     % flat couplant
        else type = '3d_b'; end                         % regular 3D
    end
    obj = clBackprojection(type);
    
    % Setup
    obj.Fs = par.fs;
    obj.RoI = par.limits;
    obj.N = par.n;
    obj.Nz = par.nZ;
    if isempty(par.limit_sensors)
        obj.Projections = par.proj;
    else
        obj.Projections = numel(par.limit_sensors);
    end
    obj.TimeResolution = par.timeres;
    obj.Sensor = {par.pos_sensor};
    obj.ImpulseResponse = par.impresp;
    obj.FrequencyVector = par.freqVec;
    obj.AmplitudeVector = par.ampVec;
    
    % Depth Correction Parameters
    if datainfo.HWDesc.AxialOffset ~= 0, obj.AxisOffset = double(datainfo.HWDesc.AxialOffset) ; end
    if datainfo.HWDesc.LightSpotSize~= 0, obj.LightSpotSize = double(datainfo.HWDesc.LightSpotSize) ; end
    
    % Couplant
    if par.CouplantSoS ~= 0, obj.CouplantSoS = par.CouplantSoS; end
    if par.CouplantRadius ~= 0, obj.CouplantRadius = par.CouplantRadius; end
    if par.CouplantOrigin ~= 0, obj.CouplantOrigin = par.CouplantOrigin; end
    
    obj.Setup;
%     obj

  
else
    error('Only backprojection (direct) supported at this time');
end


%% reconstruct
R = zeros(par.n,par.n,par.nZ,size(selMat,1),'single');
if strcmp(par.image_select,'direct'),
    sigMatF = zeros([2048 size(sigMat,2) size(sigMat,3)],'single');
end
wls = zeros(size(selMat,1),1);
zpos = zeros(size(selMat,1),1);
ts = zeros(size(selMat,1),1);

% progress bar
if (par.progress),
    wbar = waitbar(0,'Estimated Time Left: Unknown','Name','Reconstruction','CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
end

currentpath = cd('..');
parentpath = pwd();
newpath = fileparts(which('nnls_cg_armijo_reg_gpu'));
%fprintf('%s -> %s\n',currentpath,newpath);

tim = tic;
for jj = 1:size(selMat,1)
    id = selMat(jj);
    
    if (isnan(id)), continue; end;
    
    if logical(par.progress) && ~isempty(logical(getappdata(wbar,'canceling')))
        close(wbar);
        delete(wbar);
        error('User cancelled the reconstruction');
    end

    % Recon
    
    if strcmp(par.image_select,'direct') || strcmp(par.image_select,'model_lin') || strcmp(par.image_select,'nonNeg_MB'),
    
        % REMARK: mua scaling to make it fit expected 1/cm result - doublecheck
        % after bugfix in ViewMSOT
        res = obj.Run(sigMat(:,:,jj), scalefac(jj), par.c, par.ua(jj)*100);
        R(:,:,:,jj) = res;
        sigMatF(:,:,jj) = obj.RECON_RESULT.FilteredSig;
   
    end
    
    % Re-do reconstruction
    if strcmp(par.image_select,'model_lin'),
        
        currentpath = cd('..');
        parentpath = pwd();
        newpath = fileparts(which('lsqr_gpu')); 
        
       
        len = 2030; % number of samples
        t = [0:1/par.fs:(len-1)/par.fs]; % sampling instants
        n_proj = size(sigMatF,2);
        sizeT = length(t);
        b_vec = reshape(sigMatF(1:2030,:),sizeT*n_proj,1);
        tol = 1e-6;
        lsqr_iter = 20;
        %
        cd(newpath); % change path with kernel
        load dINT.mat;
        
        if par.lsqr_reg
            % regularization
            lambda = 1e3;
            nn = par.n*par.n;
            FILT = lambda*calculate_matrix_highpass_1(par.n);
            L = sparse(FILT);
            [ReconMBvec] = lsqr_gpu_reg(par.c, par.n, par.roi, par.pos_sensor(:,1), -par.pos_sensor(:,2), t, dI_mat,b_vec,lambda,L,lsqr_iter,tol);
        else
            [ReconMBvec] = lsqr_gpu(par.c, par.n, par.roi, par.pos_sensor(:,1), -par.pos_sensor(:,2), t, dI_mat,b_vec,tol,lsqr_iter);
        end
        
        R(:,:,:,jj) = reshape(ReconMBvec,[par.n par.n])*1e+8; 
    
    elseif strcmp(par.image_select,'nonNeg_MB'),
        
        % compute angle vector
        num_el = size(par.pos_sensor,1);
        [azimuth,elevation,r_sensor] = cart2sph(par.pos_sensor(:,1),par.pos_sensor(:,2),par.pos_sensor(:,3));
        pos_elements_arc_sph=zeros(num_el,3);
        pos_elements_arc_sph(:,2)=azimuth + deg2rad(180);
        pos_elements_arc_sph(:,3)=r_sensor;
        par.angle_sensor = pos_elements_arc_sph(:,2);

        % compute time vector and limits
        Dxy = par.roi/(par.n-1); % sampling distance in x and y
        if (isempty(par.limits)),
        par.limits(1) = par.center(1) - par.roi/2;              % x low
        par.limits(2) = par.center(1) + par.roi/2;        % x high
        par.limits(3) = par.center(3) - par.roiZ / 2;    % FoV start
        par.limits(4) = par.center(3) + par.roiZ / 2;    % FOV end

        par.limits(1) = par.center(1)-(par.roi+2*Dxy)*sqrt(3)/2; % limits for the signal
        par.limits(2) = par.center(1)+(par.roi+2*Dxy)*sqrt(3)/2; % limits for the signal
        end

        time_res = 3;
        dx = par.roi/par.n; % increment in x
        dt = dx/(time_res*par.c); % increment in t employed to make the model-based reconstruction
        len = 2030; % number of samples
        tss = [0:1/par.fs:(len-1)/par.fs]; % sampling instants
        dt = Dxy/(time_res*par.c); % increment in t employed to make the model-based
        fac = par.fs/par.c;
        pos_start = max(1,int32((par.limits(1))*fac));
        pos_end = min(len,int32((par.limits(2))*fac));
        t = double(tss(pos_start:pos_end));

        % regularization
        lambda = 1e3;
        nn = par.n*par.n;
        FILT = lambda*calculate_matrix_highpass_1(par.n);
        L = sparse(FILT);

     
        clear sigMatF1
        waitbar(jj/size(selMat,1),wbar,'Reconstruction...');
        sigMatF1 = sigMatF(1:2030,:,jj);
        
        % crop
        sigMat_MB = sigMatF1(pos_start:pos_end,:);

        % remove nans
        sigMat_MB(isnan(sigMat_MB))=0;

        % flatten
        sigMat_MB = reshape(sigMat_MB,[size(sigMat_MB,1)*size(sigMat_MB,2) 1]);
        b_vec = sigMat_MB(:);

        % Initialization
        c0 = single(par.c);
        n = single(par.n);
        image_width = single(par.roi);
        x_sensor = single(par.pos_sensor(:,1));
        y_sensor = single(par.pos_sensor(:,2));
        t = single(t);
        cd(newpath);
        load dINT.mat;
        dI_mat = single(dI_mat);
        b_vec= single(b_vec); 
        zetaRel=0.01;
        nIter=5;
        nInner=3;

        ReconMBNN0 = nnls_cg_armijo_reg_gpu(c0, n, image_width,x_sensor, y_sensor, t, dI_mat,b_vec,lambda,L,single(zeros(n*n,1)),zetaRel,nIter,nInner);
        res = reshape(ReconMBNN0,par.n,par.n);
        res = flipud(res);
        R(:,:,:,jj) = res*1e+8;

    end
    
    % time estimation
    if (par.progress),
        perimg = toc(tim)/jj;
        est = round(perimg*(size(selMat,1)-jj)); est_unit = 's';
        if (est > 120) est = round(est / 60); est_unit = 'min'; end;    
        wbar = waitbar(jj/size(selMat,1),wbar,sprintf( 'Estimated Time Left: %i%s', est, est_unit )) ;
    end
    
    % Meta information
    if (~isempty(datainfo)),
        wls(jj) = datainfo.ScanFrames(id).Wavelength;
        zpos(jj) = datainfo.ScanFrames(id).ZPos;
        ts(jj) = datainfo.ScanFrames(id).RelTime;
    end
end
toc(tim)
cd(currentpath);
if (par.progress), close(wbar);delete(wbar); end

 if strcmp(par.image_select,'direct'),
delete(obj);
clear obj;
 end
 
setenv('ITHERAFILES1',oldcldir);

%% Mirroring
R = R(:,par.n:-1:1,:,:); % mirror per default
if (par.mirror),
    if (islogical(par.mirror)),
        R = R(par.n:-1:1,par.n:-1:1,:,:);
    else
        mi = double(par.mirror);
        if mod(mi,2) == 0,
            mi = mi - 2;
            R = R(par.n:-1:1,:,:,:);
        end
        if (double(par.mirror) == 1),
            R = R(:,par.n:-1:1,:,:);
        end
    end            
end

%% Final Resizing
R = reshape(R,[par.n par.n par.nZ svec]);
R = R./fac;

sigMatF = reshape(sigMatF,[size(sigMatF,1), size(sigMatF,2) svec]);
sigMatF = sigMatF ./ fac;
    
ts = reshape(ts,svec);
wls = unique(wls);
zpos = unique(zpos);
sigMat = reshape(sigMat,[size(sigMat,1) size(sigMat,2) svec ]);


%% Optional: Saving
if ~isempty(datainfo) && par.save
    savefile = [datainfo.FolderName '\' datainfo.FriendlyName ...
        '_' par.image_select...
        '_roi' num2str(par.roi*1e3,'%.1f')...
        '_n' num2str(par.n)...
        '_c' num2str(par.c)...
        '_proj' num2str(par.proj)...
        '_' num2str(round(par.filter_f(1)*1e-3),'%.0f') 'khz'...
        '-' num2str(round(par.filter_f(2)*1e-6),'%.0f') 'MHz'...
        '.mat'];
    if par.verbose, fprintf('Saving Recon to %s...\n',savefile); end;
    save(savefile,'-v7.3','R','ts','wls','zpos','datainfo','par');
end
