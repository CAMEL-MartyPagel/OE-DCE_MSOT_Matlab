classdef scanFrame < handle 
    properties
        % Link to parent
        meta@msotData;
        FrameNum@uint32;
        
        % From MSOT File in ScanFrame
        IDOffset@uint32;
        Temperature@single = single(nan);
        Wavelength@double = double(nan);
        CorrectionFactor@single = single(nan);
        ZPos@single = single(nan);
        ZPosReal@single = single(nan);
        ZNum@uint16 = uint16(nan);
        WLNum@uint16 = uint16(nan);
        ShotNum@uint16 = uint16(nan);
        Timestamp@double = double(nan);
        RelTime@double = double(0);
        DiodeReadout@single = single(nan);
        LaserEnergy@single = single(nan);
        USOffset@uint32 = uint32(nan);
        Run@uint16 = uint16(nan);
        Repetition@uint16 = uint16(nan);
        
        % Behavioural
        cache@logical = false;
        selectedDetectors@uint16 = uint16([]);
    end
    
    properties ( SetAccess = private )
        rawSIG@uint16 = uint16([]);
        sigLoaded@logical = false;
    end
       
    properties ( Dependent = true )
        SIG@single;
        daqSIG@single;          % signal in mV, not energy corrected
        binFile@char;
        numSamples@char;
        numDetectors@char;
    end
    
    methods
        % *****************************************************************
       %% CONSTRUCTOR
        function obj = scanFrame(varargin)
            obj@handle;                 % superclass constructor
            if nargin == 0, return; end
            
            % Initialiation with Parent Object and Frame Number and
            % IDOffset
            if nargin >= 3,
                obj.meta = varargin{1};
                obj.FrameNum = varargin{2};
                obj.IDOffset = varargin{3};
            end
            if nargin >= 4,
                obj.cache = varargin{4};
            end
        end
        
        % *****************************************************************
       %% Methods
       % get laser energy used for correction
       function le = correctionEnergy(frame)
            le = frame.LaserEnergy;
            ale = frame.averageEnergy;
            if le <= ale*0.15, le = ale; end
       end         
       
       function ale = averageEnergy(frame)
            if ~isvalid(frame.meta), error('MSOT:scanFrame:NoMETAClass','ScanFrame cannot be used in this way without META Class'); end
            ale = frame.meta.MeasurementDesc.FactoryLaserEnergyTable(frame.meta.Wavelengths == frame.Wavelength);
       end
        
        % ******************************************************************
       %% Get Accessors
        function sig = get.rawSIG(obj)
            if obj.cache 
                % fill cache
                if ~obj.sigLoaded, 
                    obj.rawSIG = loadSignal(obj); 
                    if ~isempty(obj.rawSIG), obj.sigLoaded = true; end;
                end
                % access cache
                sig = obj.rawSIG; 
            else
                sig = loadSignal(obj);
            end
        end
        
        function sig = get.SIG(obj)
            sig = single(obj.rawSIG) ./ obj.CorrectionFactor;
        end

        function sig = get.daqSIG(obj)
            if ~isvalid(obj.meta), error('MSOT:scanFrame:NoMETAClass','ScanFrame cannot be used in this way without META Class'); end

            le = obj.LaserEnergy;
            ale = obj.meta.MeasurementDesc.FactoryLaserEnergyTable(obj.meta.Wavelengths == obj.Wavelength);
            if le <= ale*0.15, le = ale; end
            sig = single(obj.rawSIG) ./ obj.CorrectionFactor .* le .* 4e-3;
        end
        
        function x = get.binFile(obj)
           if ~isvalid(obj.meta), error('MSOT:scanFrame:NoMETAClass','ScanFrame cannot be used in this way without META Class'); end
           x = obj.meta.FileName;
        end

        function x = get.numSamples(obj)
           if ~isvalid(obj.meta), error('MSOT:scanFrame:NoMETAClass','ScanFrame cannot be used in this way without META Class'); end
           x = obj.meta.MeasurementDesc.RecordLength;
        end
        
        function x = get.numDetectors(obj)
           if ~isvalid(obj.meta), error('MSOT:scanFrame:NoMETAClass','ScanFrame cannot be used in this way without META Class'); end
           x = obj.meta.HWDesc.NumDetectors;
        end
        
        % ******************************************************************
       %% Set Accessors
        function set.cache(obj,val)
           obj.cache = val;
           % clear cache
           if ~val, obj.rawSIG = uint16([]); end
        end
        
        
        
        % ******************************************************************
       %% Signal Loading
        function sig = loadSignal(obj)
            sig = uint16([]);
            if isempty(obj.binFile), error('MSOT:scanFrame:NoBinFile','No binary File specified'); end
            if ~exist(obj.binFile,'file'), error('MSOT:scanFrame:BinFileNotExisting','Binary File %s does not exist',obj.binFile); end
            
            FID = fopen(obj.binFile,'r');
            if FID == -1, error('MSOT:scanFrame:BinFileOpenError','Error opening Binary File (%s)',obj.binFile); end
            
            try
                fpos = double(obj.IDOffset)*double(obj.numSamples)*double(obj.numDetectors)*2;
                ndet = obj.numDetectors;
                % adjust position and number of detectors if reduced number
                % of detectors requested
                if (~isempty(obj.selectedDetectors)), 
                    fpos = fpos+double(obj.selectedDetectors(1)-1)*double(obj.numSamples)*2; 
                    ndet = diff([obj.selectedDetectors(1) obj.selectedDetectors(end)])+1;
                end
                fseek(FID,fpos,-1);
                % read
                sig = uint16(2^16-1)-uint16(fread(FID,[obj.numSamples ndet],'uint16'));
                % truncate if selected detectors are scattered
                if (max(diff(obj.selectedDetectors)) ~= 1), 
                    sig = sig(obj.selectedDetectors-obj.selectedDetectors(1)+1,:); 
                end;
            catch ex
                warning(['Cannot read Signals for ScanFrame ' num2str(obj.FrameNum) ', skipping... (%s)'],ex.message);
            end
            
            fclose(FID);
        end
    end
end