function [p_t, t, angle_sensor] = preprocessSignals(p_t0,datainfo,varargin)

par = struct();
par.filter_f = [50 7000]*1e3;
par.proj = datainfo.HWDesc.NumDetectors*2 - 1;
par.n = 200;
par.roi = 20e-3;
par.invert = 1;
par.angle_sensor = [];
par.angle_offset = 0;
par.limitSensors = [];

par.useGPU = 0;         % use GPU - will skip most preprocessing

par.f_HPFOrder = 4;
par.f_HPFRipple = 0.01;
par.f_HPFzp = 1;

% default speed of sound from temperature
T = datainfo.AverageTemperature;
par.c = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;


% Copy parameters from input struct
if numel(varargin) >= 1
    cpar = varargin{1};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

% Derived Parameters
R = datainfo.HWDesc.Radius;
if (isempty(par.angle_sensor))
    if isfield(datainfo.HWDesc,'ProjectionData'),
        angle_sensor = datainfo.HWDesc.ProjectionData;
    else
        angle_sensor = linspace(datainfo.HWDesc.StartAngle,datainfo.HWDesc.EndAngle,datainfo.HWDesc.NumDetectors)+par.angle_offset/360*pi*2;
    end
else
    angle_sensor = par.angle_sensor;
end
fs = datainfo.HWDesc.SamplingFrequency;
c = par.c;
n = par.n;
ROI = par.roi;
coreNum = feature( 'numCores' ) ;

% invert signals if requested
if (par.invert)
    p_t0 = -p_t0;
end

%% truncate date for limitSensors
if ~isempty(par.limitSensors)
    if isempty(par.angle_sensor)
        angle_sensor = angle_sensor(par.limitSensors);
    end
    p_t0 = p_t0(:,par.limitSensors,:);
end

%% skip rest for GPU processing
if (par.useGPU),
    p_t = p_t0;
    t = [];
    return
end

%% time and space vectors
t_full = 1/fs:1/fs:size(p_t0,1)/fs; % initial full time vector
ROI_max = sqrt(2)*ROI;          % diagonal width from corner to corner
t_cor = round(R/c*fs);          % sample of COR
t_min = (ROI_max/c*fs);         % minimum number of samples for ROI
x=linspace(-ROI/2,ROI/2,n);   % x vector in distance from COR

%% filter the signal
% fprintf(['Filtering %i Signals (Highpass: %i kHz, Lowpass: %.1f MHz)\n'],size(p_t0,3),par.filter_f(1)*1e-3,par.filter_f(2)*1e-6);
% design filters
if par.filter_f(2)
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * par.filter_f(2)/fs * .9 ) ;
end
if par.f_HPFOrder == 4
    [b_HPF,a_HPF] = cheby1( 4, .01, 2 * par.filter_f(1)/fs * 1.46, 'high' ) ;
else
    [b_HPF,a_HPF] = cheby1( par.f_HPFOrder, par.f_HPFRipple, 2 * par.filter_f(1)/fs , 'high' ) ;
end

cut = 150;          % cut the first 150 samples
p_t0_filt = zeros(size(p_t0));
for jj = 1:size(p_t0,3)
    % low pass
    if par.filter_f(2)
        p_t0_filt(:,:,jj) = FiltFiltM( b_LPF, a_LPF, p_t0(:,:,jj), 1, coreNum ) ;
    else
        p_t0_filt(:,:,jj) = p_t0(:,:,jj);
    end
    
    % crop first samples and remove DC
    p_t0_filt(:,:,jj) = p_t0_filt(:,:,jj) - ones( size( p_t0_filt, 1 ), 1 ) * p_t0_filt( cut+1, :, jj ) ;
    p_t0_filt( [1:cut], : , jj ) = 0 ;
    
    % high pass
    if par.f_HPFzp
        p_t0_filt(:,:,jj) = FiltFiltM( b_HPF, a_HPF, p_t0_filt(:,:,jj), 1, coreNum ) ;
    else 
        p_t0_filt(:,:,jj) = FilterM( b_HPF, a_HPF, p_t0_filt(:,:,jj)) ;
    end
end
clear jj p_t0;



%% cutting the time and signal vector
% crop time vector
t_start = round(t_cor - t_min/2*1.2);
t_end = round(t_cor + t_min/2*1.2);
t=t_full(t_start:t_end);             % crop the time vector
clear t_start t_end;


% interpolated angle vector
if size(angle_sensor,1) == 1,
    dphi =  (angle_sensor(end)-angle_sensor(1) ) / (par.proj-1) ;
    angle_sensor_int = angle_sensor(1) : dphi : angle_sensor(end) ;
    clear dphi;
else
    angle_sensor_int = angle_sensor;
    if par.proj > size(angle_sensor,1),
        warning('No sensor interpolation performed with datasets that have ProjectionData specified');
    end
end

% interpolate and crop data
[X,Y] = meshgrid(angle_sensor,t_full);
[XI,YI] = meshgrid(angle_sensor_int,t);
for jj = 1:size(p_t0_filt,3)
    p_t(:,:,jj) = interp2(X,Y,p_t0_filt(:,:,jj),XI,YI); 
%     p_t(:,jj) = reshape(p_tt,size(p_tt,1)*size(p_tt,2),1);
end
angle_sensor = angle_sensor_int;
clear X Y XI YI angle_sensor_int jj p_tt; 
clear t_min t_cor ROI_max t_full;