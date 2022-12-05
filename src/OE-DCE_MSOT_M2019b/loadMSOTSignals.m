function [sig, datainfo, selMat, par, sigf] = loadMSOTSignals(varargin)
% LOADMSOTSIGNALS  Load binary scan data from .bin 
%   [sigMat, datainfo] = LOADMSOTSIGNALS()  
%        uses Open File Dialog to choose .msot file
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename)  
%        uses filename and path .msot file to load Data
%   [sigMat, datainfo] = LOADMSOTSIGNALS(datainfo)  
%        load signals from the scan described by datainfo
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename/datainfo,selMat)  
%        restricts loading of data to specific frame numbers supplied in
%        selMat. If empty, datainfo.ScanStructure is used. Output
%        will have the same structure as selMat.
%   [sigMat, datainfo] = LOADMSOTSIGNALS(filename/datainfo,selMat,params)  
%        supplies an additional parameter struct with the following fields:
%        - average: Average data while loading (default: 1)
%              NOTE: Only applies if data was acquired without averaging
%        - usePower: Laser Energy Handling (default: 1)
%            * 0: load plain binary information (will be scaled to fit
%                 16bit). IMPORTANT: Only use this for non-averaged data!
%            * 1: undo scaling and retain laser energy correction
%            * 2: undo scaling and laser energy correction.
%                 WARNING: requires subsequent laser energy correction
%            * 3: undo scaling and laser energy correction, then recorrect
%                 using Average Laser Energy
%        - shotsel: Array of single shots to average. Empty: all (default)
%              NOTE: Only applies on non-averaged data.
%        - filter: Filter the data after loading (default: 0)
%            * 0: Do not filter
%            * 1: Filter (zero phase filtering, i.e. forward and backward)
%            * 2: Filter (non-zero phase filtering, i.e. only forward)
%            * 3: Butter filter Order 3
%            * 4: Butter filter Order 4
%        - f_LPF: Low-pass edge frequency in Hz (default: 8e6)
%        - f_HPF: High-pass edge frequency in Hz (default: 50e3)
%        - f_HPFOrder: Order of Chebychev I HP Filter (default: 4)
%        - f_HPFRipple: Allowed passband ripple at edge freq (default: 0.1)
%
%
% Return Values:
%   1: sigMat   a matrix with dimensions equivalent to
%               datainfo.ScanStructure or selMat if supplied
%
% Example:
%   params.filter = 1;          % filter signals
%   sigmat = loadMSOTSignals('Scan_1\Scan_1.msot',[],params);
% 


selMat = [];
sig = [];
datainfo = [];
par.shotsel = [];           % select shots to load, empty is all
par.usePower = 1;           % undo scaling or apply power correction
par.f_LPF = 8000000;        % filter freq. defaults
par.f_HPF = 50000;          % filter freq. defaults
par.f_HPFOrder = 4;
par.f_HPFRipple = 0.1;
par.filter = 0;             % dont filter, 1: zerophase, 2:non-zp
par.average = 1;            % but average per default
par.useSingle = false;
par.progress = true;
par.avgdim = 0;             % dimension along which to average
par.avgcorr = 0;            % correlation limit for averaging (dropping non-correlating)
par.watercorr = 0;         % do not do correction of the data for the water absorption
par.noisewin = [];          % window for noise measurement

fs = 4e7;

% get function parameters
if nargin == 0
    datainfo = loadMSOT();
end
if (nargin >= 1)
    if (isstruct(varargin{1}) || isa(varargin{1},'msotData'))
        datainfo = varargin{1};
    else
        datainfo = loadMSOT(varargin{1});
    end
end
if nargin >= 2
    selMat = varargin{2};
end
if nargin >= 3
    cpar = varargin{3};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

% check parameters
if (isempty(datainfo))
    warning('Metadata not found, exiting...');
    return;
end
% per default, load complete dataset
if (isempty(selMat)), selMat = datainfo.ScanStructure; end

%% Progress Bar
if par.progress, wbar = waitbar(0,'Opening file...'); end;

%% open File
FID = fopen(datainfo.FileName,'r');

%% determine read structure

% select shots
if ( ~isempty(par.shotsel) && ndims(selMat) >= 5 )
   selMat = selMat(:,:,:,:,par.shotsel);
end


svec = size(selMat);        % dimension vector
% if not averaged, read structure as is in 1D
if ~par.average
    frameids = reshape(selMat,prod(svec),1);
    numavg = 1;
% otherwise put ids to average in 2D array
else
    % use last (5th) number for number of averages and shorten dim vector
    if numel(svec) == 5   
        numavg = svec(5); 
        svec = svec(1:4);       % cut size vector for output
    % otherwise the dataset is already averaged
    else, numavg = 1; 
    end;
    
    frameids = reshape(selMat,prod(svec),numavg);
    selMat = selMat(:,:,:,:,1);
end

%% data type
dtype = 'uint16';                           % default to native uint16
if numavg > 15,  dtype = 'uint32'; end       % for more averages use uint32
% if filtering or power correction required, double is necessary
if (par.filter || par.usePower), dtype = 'double'; end
if par.useSingle, dtype = 'single'; end;
    
%% initialize filters
if (par.filter && par.filter <= 2)
    if par.f_LPF, [b_LPF,a_LPF] = cheby1( 8, .01, 2 * par.f_LPF/fs * .9 ) ; end
    
    if par.f_HPF 
        if par.f_HPFOrder == 4
            [b_HPF,a_HPF] = cheby1( 4, .01, 2 * par.f_HPF/fs * 1.46, 'high' ) ; 
        elseif par.f_HPFOrder == 1
            [b_HPF,a_HPF] = cheby1( 1, 1, 2 * par.f_HPF/fs , 'high' ) ; 
        end
    end
end
% Alex Special Filter for PreFiltering
if (par.filter && par.filter == 3)
    n=4; 
    Wn=[par.f_HPF/(fs/2) par.f_LPF/(fs/2)]; 
    ftype='bandpass';
    [b_LPF,a_LPF] = butter(n,Wn,ftype);
end

if (par.filter && par.filter == 4)
   [b_LPF,a_LPF] = butter(3,2*[par.f_HPF par.f_LPF]/fs);
end

%% read data
framenum = size(frameids,1)-1;
if par.progress, waitbar(0,wbar,'Loading Signals...'); end
numproj = datainfo.HWDesc.NumDetectors;
numsamples = double(datainfo.MeasurementDesc.RecordLength);
sig = zeros(numsamples,numproj,size(frameids,1),dtype);
if ~isempty(par.noisewin) || par.avgdim > 0 
   sig_raw = zeros(numsamples,numproj,size(frameids,1),'uint16');
end
        
for j = 1:size(frameids,1)
    if par.progress, waitbar(j/framenum,wbar); end
    tmp = zeros(numsamples,numproj,dtype);
    readcount = 0;          % counter for read frames
    % iterate through second dimension to capture averages
    for ii = 1:size(frameids,2)
        % get current id and frame data
        id = frameids(j,ii);
        if (isnan(id)), continue; end;
        frame = datainfo.ScanFrames(id);
        offset = frame.IDOffset;
%         fprintf('id: %i - offset: %i\n',id,offset);
        
        readcount = readcount + 1;

        % position and read uint16 data
        tmp2 = zeros(numsamples,numproj,dtype);
        try
            fseek(FID,double(offset)*numsamples*numproj*2,-1);
            tmp2(:,:) = fread(FID,[numsamples numproj],'uint16');
            tmp2(1:2,:) = 0;
        catch 
            warning(['Cannot Read ScanFrame ' num2str(id) ', skipping...']);
            readcount = readcount - 1;
        end

        % Correction Factors
        try
            le = frame.LaserEnergy;
            ale = datainfo.MeasurementDesc.FactoryLaserEnergyTable(datainfo.Wavelengths == frame.Wavelength);
            if le <= ale*0.15, le = ale; end

            cf = frame.CorrectionFactor;
            if cf == 0
                warning('Correction Factor is 0 for frame %i - this could indicate a malfunction. Please observe data scaling carefully',ii);
                cf = 1;
                if frame.LaserEnergy == 0  % in this case the calibration was not there, so no laser energy correction was done, so we should devide by average energy 
                    if ale == 0,  % bugfix for the case of absent calibration value - better than correcting for nothing
                        ale = mean(datainfo.MeasurementDesc.FactoryLaserEnergyTable);
                        warning('No Laser Energy correction for wavelength %i - this could indicate a malfunction. Please observe data scaling carefully, now correcting with the average energy of all wavelengths (%.1fmJ)',frame.Wavelength, ale);                        
                    end
                    cf = ale;
                end
            end
        catch 
            error('Missing or faulty META-DATA, please contact iThera Support');
        end
        
        % correct for scaling
        try 
            tmp3 = tmp2 ./ cf .* le;
        catch
            tmp3 = tmp2;
        end
        
        if par.usePower
            if (par.usePower >= 2)
                tmp2 = tmp2 ./ cf .* le;
                % Re-Correct with Average Laser Energy
                if par.usePower == 3
                    tmp2 = tmp2 ./ ale;
                end
            elseif (par.usePower <= -1) && ~isempty(par.noisewin)
                % bring values back to native DAQ scaling, ignoring
                % everything the SW did
                ntmp = tmp2(par.noisewin,:);
                tmp2 = cast(double(tmp2)./mean(ntmp(:)).*2048,dtype);
                % normalise per channel instead
                %tmp2 = cast(double(s)./(ones(2030,1)*mean(s(noisewin,:))).*32768,getfield(whos('s'),'class'));
                
                % do laser energy correction
                if par.usePower == -1          % per pulse
                    tmp2 = tmp2 ./ le;
                elseif par.usePower == -3     % average
                    tmp2 = tmp2 ./ ale;
                end
            else
                tmp2 = tmp2 ./ cf;
            end
        end
        
        
        % filtering (1,3,..: filter, 2,4,...: zero phase)
%         if (par.filter) && (mod(par.filter,2) == 1)
        if (par.filter) && (mod(par.filter,2) == 1)
            if par.f_LPF, tmp2 = FiltFiltM( b_LPF, a_LPF, tmp2, 1, 2 ) ; end
            if (par.filter == 1) && par.f_HPF, tmp2 = FiltFiltM( b_HPF, a_HPF, tmp2, 1, 2 ) ; end
        elseif par.filter == 2
            if par.f_LPF, tmp2 = FilterM( b_LPF, a_LPF, tmp2 ) ; end
            if par.f_HPF, tmp2 = FilterM( b_HPF, a_HPF, tmp2 ) ; end
        elseif par.filter == 4
            tmp2 = filtfilt(b_LPF,a_LPF,double(tmp2));
        end
        
        % add to vector
        tmp = tmp + tmp2;
        
    end
    % average all if any frames were read (otherwise leave zeros)
    if readcount
        sig(:,:,j) = tmp ./ readcount;
        if ~isempty(par.noisewin) || par.avgdim > 0 
            sig_raw(:,:,j) = tmp3 ./ readcount;
        end
    end
end;
sig = reshape(sig,[numsamples numproj svec]);
if ~isempty(par.noisewin) || par.avgdim > 0 
    sig_raw = reshape(sig_raw,[numsamples numproj svec]);
end
clear ii j tmp tmp2 tmp3 frame id;

%% Correction of the data for the water absorption
if par.watercorr == 1
    WaterAbsorptionCoeff = datainfo.MeasurementDesc.WaterAbsorptionCoeff;
    PathLengthInWater = datainfo.MeasurementDesc.PathLengthInWater;
    T = exp( - WaterAbsorptionCoeff * (PathLengthInWater)) ;
    sig_tmp = reshape(sig,[size(sig,1) size(sig,2) size(sig,3)*size(sig,4)*size(sig,5) size(sig,6)]);

    for t = 1:size(sig_tmp,3)
        for wl = 1:size(sig_tmp,4)
            sig_tmp(:,:,t,wl) = sig_tmp(:,:,t,wl)/T(wl,1) ;
        end
    end
    sig = reshape(sig_tmp,[numsamples numproj svec]);
end

%% average if required
if par.avgdim > 0 
    avgsel = 1:size(sig,par.avgdim+2);  % per default, take all frames
    if par.avgcorr > 0 
%         par.avgcorr
        % ** Generate Mean Filtered Signal
        [b_HPF,a_HPF] = cheby1( 1, 1, 2 * 500e3/fs , 'high' ) ;
        [b_LPF,a_LPF] = cheby1( 8, .01, 2 * 8e6/fs * .9 ) ;
%         sigf = squeeze(sig);
        if exist('sig_raw','var') && nnz(isnan(sig_raw(:))) == 0
            sigf = double(reshape(sig_raw,[size(sig,1) size(sig,2) size(sig,3)*size(sig,4)*size(sig,5) size(sig,6)]));
        else
            sigf = double(reshape(sig,[size(sig,1) size(sig,2) size(sig,3)*size(sig,4)*size(sig,5) size(sig,6)]));
        end            
        if numel(size(sigf)) > 4, error('too many dimensions for selective averaging'); end
        for t = 1:size(sigf,3)
            for wl = 1:size(sigf,4)
                sigf(:,:,t,wl) = FilterM( b_LPF, a_LPF, sigf(:,:,t,wl) ) ;
                sigf(:,:,t,wl) = FilterM( b_HPF, a_HPF, sigf(:,:,t,wl) ) ;
            end
        end
        sigwindow = 600:1500;

        % in case there is something to average
        if size(sig,par.avgdim+2) > 1
            sigm = squeeze(mean(sigf,4));
            sigm = (reshape(sigm(sigwindow,:,:),[numel(sigwindow)*size(sigm,2) size(sigm,3)]));
            R = corrcoef(sigm);

%             cref = round(median(avgsel));
            corr = R(:,end);
            avgsel = avgsel(corr > par.avgcorr);
            par.avgcorr = mean(corr(corr >= par.avgcorr));
        else
            par.avgcorr = 1;
        end
    end
    par.avgretain = numel(avgsel);  

    % average frames
    sig = mean(sig(:,:,:,:,avgsel,:,:),par.avgdim+2);
    selMat = min(selMat,[],par.avgdim);
%     datainfo.LaserEnergy = mean(datainfo.LaserEnergy(:,:,avgsel,:),par.avgdim);
    
    % calculate correlation of averaged frames in the WL dimension
    if par.avgcorr > 0 && par.avgdim == 3
      sigf = mean(sigf(:,:,avgsel,:),3);
      if size(sigf,4) > 1
        CC = corrcoef(reshape(squeeze(sigf(sigwindow,:,:,:)),[numel(sigwindow)*size(sigf,2) size(sigf,4)]));
        par.avgcons = mean(CC(:,round(size(sigf,4)/2+1)));
      else
        par.avgcons = 1;
      end
    end
       
end


%% noise calculation
% TODO: Actually, this should always be done on unscaled signals in V)
if exist('sigf','var') && ~isempty(par.noisewin)
    sn = shiftdim(std(sigf(par.noisewin,:,:,:),0,1),1);
    par.noise = sn.*4;         % STD in µV 
end

% hide dialog
if par.progress, close(wbar); end

%% close file
fclose(FID);