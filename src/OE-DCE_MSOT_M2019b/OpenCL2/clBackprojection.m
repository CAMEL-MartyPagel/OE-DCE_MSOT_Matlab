% MATLAB class wrapper to the underlying C++/openCL functionality
classdef clBackprojection < handle
    properties (SetAccess = private, Hidden = true)
        % Handle to the underlying C++ class instance
        m_objectHandle; 
        % string that defines the reconstruction type:
        % '2d'     2d backprojection
        % '3d_b'   3d backprojection no couplant
        % '3d_f'   3d backprojection with flat couplant
        % '3d_c'   3d backprojection with curved couplant
        m_mode;
    end
    properties (SetAccess = private)
        SETUP_PARAMS = struct(...
            'ImpulseResponse', [], ...
            'SensitivityMap', [], ... % not implemented yet
            'SensorRadius', [], ...
            'SensorAngle', [], ...
            'Filter', [0,0], ...
            'Fs', 4e7, ...
            'SoS', 0, ...
            'CouplantSoS', 0, ...
            'CouplantRadius', 0, ...
            'CouplantOrigin', 0, ...
            'TimeResolution', 3, ...
            'ROILow', 0,... 
            'ROIHigh', 0,...
            'ROIZLow', 0,...
            'ROIZHigh', 0,...
            'RotationAngle', 0, ...
            'AThresh', 0, ... % for wavelets
            'Thresh', 0, ...  % for wavelets
            'AxisOffset', 0, ...
            'LightSpotSize', 0, ...
            'Channels', 0, ...
            'Projections', 0, ...
            'Samples', 2030, ...
            'InterpolationType', 3, ...
            'LsqrIter', 20, ...
            'N', 0, ...
            'Nz', 0, ...
            'FilterType', [0,0], ...
            'Zerophase', [0,0], ...
            'Order', [0,0], ...
            'rp', [0,0], ...
            'rs', [0,0], ...
            'FrequencyMultiplier', [0,0], ...
            'FrequencyVector', [], ...
            'AmplitudeVector', []);

        RECON_PARAMS = struct(...
            'PowerFactor', 1, ...
            'SoS', 0, ...
            'ua', 0, ...
            'XY', [], ...
            'YZ', [], ...
            'ZX', [], ...
            'FilteredSig', [], ...
            'LimitsX', [0,0], ...
            'LimitsY', [0,0], ...
            'LimitsZ', [0,0], ...
            'XYIntensityRange', [0,0], ...
            'YZIntensityRange', [0,0], ...
            'ZXIntensityRange', [0,0]);
        RECON_RESULT = struct(...
            'Recon', [], ...
            'FilteredSig', [], ...
            'XY', [], ...
            'YZ', [], ...
            'ZX', [], ...
            'XYIntensityRange', [0,0], ...
            'YZIntensityRange', [0,0], ...
            'ZXIntensityRange', [0,0]);
    end
    properties (Dependent = true)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SETUP PARAMS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % impulse response vector
        ImpulseResponse;
        % sensor description
        % a cell array containing either an array of cartesian positions
        % or 
        % a radius and a vector of angles
        Sensor;
        % acquisition sampling frequency: always 40 MHz.
        Fs;
        % couplant speed of sound
        CouplantSoS;
        % couplant radius
        CouplantRadius;
        % couplant origin
        CouplantOrigin;
        % time resolution
        TimeResolution;
        % range of interest description
        RoI;
        % number of projections
        Projections;
        % Axial offset
        AxisOffset;
        % lazer light spot size
        LightSpotSize;

        % Horizontal and Vertical pixel/voxel resolution
        N;
        % Depth voxel resolution. Must be 0 in 2D.
        Nz;
        % interpolation type:
        % 1 = neaarest neighbour
        % 2 = linear
        % 3 = bicubic 
        % Always use bicubic. It performs the same with better quility.
        InterpolationType;
        
        %%%%%%%%%%%%%%%%%
        % ANALOG FILTERS 
        %%%%%%%%%%%%%%%%%
        % define low frequency and high frequensy edges of the spectral
        % response of an analog filter (if used).
        Filter;
        % defines the analog lowpass and highpass filters if used.
        % 0 none
        % 1 Butterworth
        % 2 Xhebyshev I
        % 3 Chebyshev II
        FilterType;
        % array (for lowpass and highpass), 0 = no, 1 = yes
        Zerophase;
        % array (for lowpass and highpass)
        Order;
        % array (for lowpass and highpass), pass band ripple for Chebyshev
        rp;
        %%%%%%%%%%%%%%%%%
        % DIGITAL FILTER
        %%%%%%%%%%%%%%%%%
        % If FrequencyVector is empty the above parameters are used to setup
        % preprocessing. In the case of 3D cup the freqVec is used 
        % along with ampVec. freqVeq is a vector of frequencies. 
        % The first is 0 and the last is 1. The rest of the values 
        % are actual frequencies normalized at the Nyquist frequency.
        FrequencyVector;
        % array of complex values corresponding to spectral amplitude at 
        % the frequencies of freqVec. The actual spectrum is a linear 
        % interpolation of these values at he freqVec frequencies. 
        % Usually a band limited differentiator is described in this way
        % and a FIR filter is designed.
        AmplitudeVector;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RECON PARAMS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % a cell array of 3 arrays
        % each contain the min and max coordinate to MIP on the X,Y,Z
        Limits;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OUTPUT PROPERTIES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculated 2d sensor deviation angle from vertical symmetry
        RotationAngle; 
        % reconstructed data
        Recon;
        % the filtered signal after preproseccing
        FilteredSig;
        % for 3d reconstruction
        % xy plane projection of the volumetric data
        XY;
        % yz plane projection of the volumetric data
        YZ;
        % zx plane projection of the volumetric data
        ZX;
        % a cell array of 3 arrays
        % each containing the min and max intensities of XY,YZ,ZX
        IntensityRange;
    end
    methods
        % Constructor - Create a new C++ class instance 
        % mode can be:
        % '2d'     2d backprojection
        % '3d_b'   3d backprojection no couplant
        % '3d_f'   3d backprojection with flat couplant
        % '3d_c'   3d backprojection with curved couplant
        function this = clBackprojection(mode)
            this.m_mode = mode;
            this.m_objectHandle = clBackprojectionSync(mode);
        end
        
        % Destructor - Destroy the C++ class instance
        function delete(this)
            clBackprojectionSync('d', this.m_objectHandle);
        end

        % setup internal state according to the parameters supplied
        % already in SETUP_PARAMS
        % returns RotationAngle.
        function rAngle = Setup(this)
            this.SETUP_PARAMS.RotationAngle = clBackprojectionSync('s', this.m_objectHandle, this.SETUP_PARAMS);
            rAngle = this.SETUP_PARAMS.RotationAngle;
        end

        % do the backprojection of the data according to the last setup
        % sigMat_u16    is the acquired in 16bit unsigned format.
        % PowerFactor   is the signal amplitude correction factor
        % c             is the speed of sound
        % ua            is the optional absorption coefficient
        function r = Run(this, sigMat_u16, PowerFactor, c, ua)
            if nargin == 4, ua = 0; end;
            this.RECON_PARAMS.PowerFactor = PowerFactor;
            this.RECON_PARAMS.SoS = c;
            this.RECON_PARAMS.ua = ua;
            this.RECON_RESULT = clBackprojectionSync('r', this.m_objectHandle, sigMat_u16, this.RECON_PARAMS);
            r = this.RECON_RESULT.Recon;
        end
        
        % Perform automatic speed of sound selection based 
        % on the assessment of image quality. This is a pure c++/openCL
        % implementation
        function sos = Autofocus(this, sigMat, pcf, ua, lb, ub, nFineRange, fineStep)
            sos = clBackprojectionSync('a', this.m_objectHandle, sigMat, pcf, ua, lb, ub, nFineRange, fineStep);
        end
        % Perform automatic speed of sound selection based 
        % on the assessment of image quality
        function sos = AutofocusScript(this, sigMat, pcf, ua, lb, ub, nFineRange, fineStep)
            % Author: Elena Nasonova (2013)

            par = this.RECON_PARAMS;
            par.PowerFactor = pcf;
            par.sigMat = sigMat;
            par.ua = ua;

            lb = floor(lb);
            ub = ceil(ub);

            SoS = lb:fineStep:ub;
            metrics=zeros(1,numel(SoS));

            for i=1:4:numel(SoS)
                metrics(i) = this.myCostFunction(SoS(i),par);
            end
            [~,i] = min(metrics);
            m = i+1;
            for k=1:nFineRange
                if m <= numel(SoS) && metrics(m) == 0,
                    metrics(m) = this.myCostFunction(SoS(m),par);
                end
                m = m+1;
            end
            m = i-1;
            for k=1:nFineRange
                if (m >= 1 && m <= numel(SoS)) && metrics(m) == 0,
                    metrics(m) = this.myCostFunction(SoS(m),par);
                end
                m = m-1;
            end

            metrics = sgolayfilt(metrics,1,3);

            [~,i] = min(metrics);
            sos = SoS(i);
            fprintf('SoS = %g, Bounds[%g, %g]\n', sos, lb, ub);
        end        
        
        %
        % Methods related to setup params
        %
        % range of interest description
        function set.RoI(obj,v)
            obj.SETUP_PARAMS.ROILow = v(1);
            obj.SETUP_PARAMS.ROIHigh = v(2);
            if length(v) > 2,
                obj.SETUP_PARAMS.ROIZLow = v(3);
                obj.SETUP_PARAMS.ROIZHigh = v(4);
            end
        end
        function v = get.RoI(obj)
            if isequal(obj.m_mode, '2d'),
                v = [obj.SETUP_PARAMS.ROILow, ...
                    obj.SETUP_PARAMS.ROIHigh];
            else
                v = [obj.SETUP_PARAMS.ROILow, ...
                    obj.SETUP_PARAMS.ROIHigh, ...
                    obj.SETUP_PARAMS.ROIZLow, ...
                obj.SETUP_PARAMS.ROIZHigh];
            end
        end
        %-----------------------------------------------------
        % couplant speed of sound
        function set.CouplantSoS(obj,v)
            obj.SETUP_PARAMS.CouplantSoS = v;
        end
        function v = get.CouplantSoS(obj)
            v = obj.SETUP_PARAMS.CouplantSoS;
        end
        %-----------------------------------------------------
        % couplant radius
        function set.CouplantRadius(obj,v)
            obj.SETUP_PARAMS.CouplantRadius = v;
        end
        function v = get.CouplantRadius(obj)
            v = obj.SETUP_PARAMS.CouplantRadius;
        end
        %-----------------------------------------------------
        % couplant origin
        function set.CouplantOrigin(obj,v)
            obj.SETUP_PARAMS.CouplantOrigin = v;
        end
        function v = get.CouplantOrigin(obj)
            v = obj.SETUP_PARAMS.CouplantOrigin;
        end
        %-----------------------------------------------------
        function set.AxisOffset(obj,v)
            obj.SETUP_PARAMS.AxisOffset = v;
        end
        function v = get.AxisOffset(obj)
            v = obj.SETUP_PARAMS.AxisOffset;
        end
        %-----------------------------------------------------
        function set.LightSpotSize(obj,v)
            obj.SETUP_PARAMS.LightSpotSize = v;
        end
        function v = get.LightSpotSize(obj)
            v = obj.SETUP_PARAMS.LightSpotSize;
        end
        %-----------------------------------------------------
        % time resolution
        function set.TimeResolution(obj,v)
            obj.SETUP_PARAMS.TimeResolution = v;
        end
        function v = get.TimeResolution(obj)
            v = obj.SETUP_PARAMS.TimeResolution;
        end
        %-----------------------------------------------------
        % RotationAngle is the detected deviation from the vertical 
        % of the symmetry axis of a 2d sensor
        function v = get.RotationAngle(obj)
            v = obj.RECON_RESULT.RotationAngle;
        end
        %-----------------------------------------------------
        % number of projections
        function set.Projections(obj,v)
            obj.SETUP_PARAMS.Projections = v;
        end
        function v = get.Projections(obj)
            v = obj.SETUP_PARAMS.Projections;
        end
        %-----------------------------------------------------
        % vector of frequencies.
        % The first is 0 and the last is 1. The rest of the values 
        % are actual frequencies normalized at the Nyquist frequency.
        function set.FrequencyVector(obj,v)
            obj.SETUP_PARAMS.FrequencyVector = v;
        end
        function v = get.FrequencyVector(obj)
            v = obj.SETUP_PARAMS.FrequencyVector;
        end
        %-----------------------------------------------------
        % array of complex values corresponding to spectral amplitude at 
        % the above frequencies. The actual spectrum is a linear 
        % interpolation of these values. 
        function set.AmplitudeVector(obj,v)
            obj.SETUP_PARAMS.AmplitudeVector = v;
        end
        function v = get.AmplitudeVector(obj)
            v = obj.SETUP_PARAMS.AmplitudeVector;
        end
        %-----------------------------------------------------
        % impulse response vector
        function set.ImpulseResponse(obj,v)
            obj.SETUP_PARAMS.ImpulseResponse = v;
        end
        function v = get.ImpulseResponse(obj)
            v = obj.SETUP_PARAMS.ImpulseResponse;
        end
        %-----------------------------------------------------
        % define low frequency and high frequensy edges of the spectral
        % response of an analog filter (if used).
        function set.Filter(obj,v)
            obj.SETUP_PARAMS.Filter = v;
        end
        function v = get.Filter(obj)
            v = obj.SETUP_PARAMS.Filter;
        end
        %-----------------------------------------------------
        % defines the analog lowpass and highpass filters if used.
        % 0 none
        % 1 Butterworth
        % 2 Xhebyshev I
        % 3 Chebyshev II
        function set.FilterType(obj,v)
            obj.SETUP_PARAMS.FilterType = v;
        end
        function v = get.FilterType(obj)
            v = obj.SETUP_PARAMS.FilterType;
        end
        %-----------------------------------------------------
        % array (for lowpass and highpass), 0 = no, 1 = yes
        function set.Zerophase(obj,v)
            obj.SETUP_PARAMS.Zerophase = v;
        end
        function v = get.Zerophase(obj)
            v = obj.SETUP_PARAMS.Zerophase;
        end
        %-----------------------------------------------------
        % analog filter order array (for lowpass and highpass)
        function set.Order(obj,v)
            obj.SETUP_PARAMS.Order = v;
        end
        function v = get.Order(obj)
            v = obj.SETUP_PARAMS.Order;
        end
        %-----------------------------------------------------
        % array (for lowpass and highpass), pass band ripple for Chebyshev
        function set.rp(obj,v)
            obj.SETUP_PARAMS.rp = v;
        end
        function v = get.rp(obj)
            v = obj.SETUP_PARAMS.rp;
        end
        %-----------------------------------------------------
        % sets the Horizontal and Vertical pixel/voxel resolution
        function set.N(obj,v)
            obj.SETUP_PARAMS.N = v;
        end
        % gets the Horizontal and Vertical pixel/voxel resolution
        function v = get.N(obj)
            v = obj.SETUP_PARAMS.N;
        end
        %-----------------------------------------------------
        % sets the Depth voxel resolution
        function set.Nz(obj,v)
            obj.SETUP_PARAMS.Nz = v;
        end
        % gets the Depth voxel resolution
        function v = get.Nz(obj)
            v = obj.SETUP_PARAMS.Nz;
        end
        %-----------------------------------------------------
        % sets the acquisition sampling frequency
        function set.Fs(obj,v)
            obj.SETUP_PARAMS.Fs = v;
        end
        % gets the acquisition sampling frequency
        function v = get.Fs(obj)
            v = obj.SETUP_PARAMS.Fs;
        end
        %-----------------------------------------------------
        % sensor description
        function set.Sensor(obj,v)
            if ~iscell(v),
                error ('Sensor must be a 1 or 2 cell array');
            elseif numel(v) == 1,
                obj.SETUP_PARAMS.SensorRadius = double(v{1});
                obj.SETUP_PARAMS.SensorAngle = [];
                obj.SETUP_PARAMS.Channels = numel(v{1})/3;
            elseif numel(v) == 2,
                obj.SETUP_PARAMS.SensorRadius = double(v{1});
                obj.SETUP_PARAMS.SensorAngle = double(v{2});
                obj.SETUP_PARAMS.Channels = numel(v{2});
            else
                error ('Sensor must be a 1 or 2 cell array');
            end
        end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Methods related to reconstruction params
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.Limits(obj,v)
            obj.RECON_PARAMS.LimitsX = v(1);
            obj.RECON_PARAMS.LimitsY = v(2);
            obj.RECON_PARAMS.LimitsZ = v(3);
        end
        function v = get.Limits(obj)
            v = {obj.RECON_PARAMS.LimitsX,...
                obj.RECON_PARAMS.LimitsY,...
                obj.RECON_PARAMS.LimitsZ };
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Methods related to reconstruction result
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %-----------------------------------------------------
        % reconstruction result
        function v = get.Recon(obj)
            v = obj.RECON_RESULT.Recon;
        end
        %-----------------------------------------------------
        % return a cell array of 3 arrays
        % each containing the min and max intensities of XY,YZ,ZX
        function v = get.IntensityRange(obj)
            v = {obj.RECON_RESULT.XYIntensityRange,...
                obj.RECON_RESULT.YZIntensityRange,...
                obj.RECON_RESULT.ZXIntensityRange};
        end
        %-----------------------------------------------------
        % planar Minimum intensity projections for the 3D case
        function v = get.XY(obj)
            v = obj.RECON_RESULT.XY;
        end
        function v = get.YZ(obj)
            v = obj.RECON_RESULT.YZ;
        end
        function v = get.ZX(obj)
            v = obj.RECON_RESULT.ZX;
        end
        %-----------------------------------------------------
        % the filtered signal before reconstruction
        function v = get.FilteredSig(obj)
            v = obj.RECON_RESULT.FilteredSig;
        end
        
        %-----------------------------------------------------
        % this method dows 2D reconstruction as backproject_gpu does.
        % not needed but the competition is heavy
        function [Recon,filteredSignal] = backproject(obj, sigMat,n,res,c,pos_elements,fmin,fmax,ampmin,ampmax,IRF)


            ua=0;
            sigMat_u16 = uint16(sigMat);

            % Other Parameters 
            fs = 4e7;
            pwrCorr = 1.0;
            num_channels=size(pos_elements,1);
            proj= 2*num_channels-1;

            %Caculating the limits
            image_width = res*n*1e-3;
            center = 1015*c/fs;
            limits(1) = center - image_width/2; % LimitLow
            limits(2) = center + image_width/2; % LimitHigh


            f1 = fmin;
            f2 = fmax;
            freqVec = [0, f1*0.8, f1,     f2,    f2*1.2, 20]/20;   % frequency
            ampVec =  [0, 0,         ampmin,  ampmax, 0          0]; % differentiator response

            try
                % Setup
                obj.RoI = limits;
                obj.N = n;
                obj.Nz = 0;
                obj.Projections = proj;
                obj.TimeResolution = 3;
                obj.Sensor = {pos_elements};
                obj.ImpulseResponse = IRF;
                obj.FrequencyVector = freqVec;
                obj.AmplitudeVector = ampVec;

                obj.Setup;

                % Do Recon
                Recon = obj.Run(sigMat_u16, pwrCorr, c, ua);
                filteredSignal = obj.FilteredSig;

            catch ex
                rethrow(ex) ;
            end
        end
              
    end
    methods (Access = private)
        function metric_value = myCostFunction(obj, c, par)
            % Cost function to be minimized 
            % Computes the metric related to the image focusing defined by 'type_corr'
            % Author: Elena Nasonova (2013)

            % construct current speed of sound
            par.SoS = c;
            ReconResult = clBackprojectionSync( 'r', obj.m_objectHandle, par.sigMat, par);
            R = ReconResult.Recon;
            if length(size(R)) > 2,
                R = max(R,[],3);
            end

            Gx = conv2(R, [-1, 0, 1; -2, 0, 2; -1, 0, 1], 'same');
            Gy = conv2(R, [1, 2, 1; 0, 0, 0; -1, -2, -1], 'same');
            G2 = (Gx.^2+Gy.^2);
            % calculate the focus metric      
            metric_value = -var(G2(:))/1500;
        end
  
    end    
end
