classdef msotData < handle
    properties
        verbose logical = false;
        progress logical = true; % --> Triggers use of waitbar while creating the object
    end
    
    properties ( SetAccess = private )
        DOM
        XMLFileName char; 
        FileName char;
        RealPath char;
        StudyName char;
        StudyPath char;
        
        Name char;
        Comment char;
        FriendlyName char;
        FolderName char;
        Scientist char;
        CreationTimeTxt char;
        CreationTime double {mustBeScalarOrEmpty};
        Complete logical {mustBeScalarOrEmpty} = true;        % legacy
        is3D logical {mustBeScalarOrEmpty} = false;
        is2D logical {mustBeScalarOrEmpty} = false;
        TransducerType char;
        
        MeasurementDesc struct {mustBeScalarOrEmpty} = struct(); % Does not appear when calling obj in constructor
        HWDesc struct {mustBeScalarOrEmpty} = struct();          % Does not appear when calling obj in constructor
        USpresent logical {mustBeScalarOrEmpty} = false;
        US struct {mustBeScalarOrEmpty};
        OAMPreset struct {mustBeScalarOrEmpty} = struct();
        RepNum uint16 {mustBeScalarOrEmpty};                    % Does not appear when calling obj in constructor
        RunNum uint16 {mustBeScalarOrEmpty};                    % Does not appear when calling obj in constructor
        ShotNum uint16 {mustBeScalarOrEmpty};                   % Does not appear when calling obj in constructor
        AverageTemperature single {mustBeScalarOrEmpty};        % Does not appear when calling obj in constructor
        
        ScanFrames scanFrame ;             % Does not appear when calling obj in constructor
        ScanStructure single;                    % Does not appear when calling obj in constructor
        RelTime double;                          % Does not appear when calling obj in constructor
        LaserEnergy single;                      % Does not appear when calling obj in constructor
        DiodeReadout single;                     % Does not appear when calling obj in constructor
        CorrectionFactor single;                 % Does not appear when calling obj in constructor
        Temperature single;                      % Does not appear when calling obj in constructor
        Wavelengths double  ;               
        ZPositions single  ;                % Does not appear when calling obj in constructor
        Timestamps double   ;                % Does not appear when calling obj in constructor
        idoffset uint32;        % lookup matrix
        
        ImgFileName char;
        ImgSession struct  ;
        
        ReconNode reconNode  ;
        MSPNode mspNode  ;
        
        structureTime single = single(nan);
    end

    properties ( Access = private )
        structureLoaded logical = false;         % Does not appear when calling obj in constructor --> Private
    end
    
    properties ( Dependent = true )
        RunString cell;                          % Does not appear when calling obj in constructor   
        ZPositionString cell;                    % Does not appear when calling obj in constructor
        RepString cell;                          % Does not appear when calling obj in constructor
        WavelengthString cell;
    end
    
    methods
        % *****************************************************************
       %% CONSTRUCTOR
        function obj = msotData(fileName,progress)
            % for array initialisation
            if nargin == 0, return; end
            if nargin >= 2, obj.progress = progress; end
            
            % check if file exists
            if ~exist(fileName,'file')
                error('File %s does not exist - please specify .msot file',fileName);
            end
            
            % check file ending
            tok = regexpi(fileName,'(.+)\.msot(.+)?$','tokens');
            if isempty(tok)
                error('Please specify a .msot file (%s)', fileName);
            end
            binFile = [ tok{1}{1} '.bin' ];
            
            % find directory and study name
            sepind = strfind(fileName, filesep);
            if isempty(sepind)
                dirName = [ '.' filesep ];
                studyPath = ['..' filesep];
                studyName = '<unknown>';
            else
                dirName = fileName(1:sepind(end));
                studyPath = fileName(1:sepind(end-1)-1);
                if numel(sepind) < 3
                    ss = 0;
                else
                    ss = sepind(end-2);
                end
                studyName = fileName(ss+1:sepind(end-1)-1);
            end
            
            
            % read XML
            if obj.progress, wbar = waitbar(0,'Parsing MSOT File'); end
            try
                patdoc = javaMethod( 'parse', 'com.itheramedical.msotbeans.DataModelMsotProjectDocument$Factory', java.io.File( fileName ) ) ; % java method, java class
            catch ex
                error( ['Selected MSOT file could not be loaded. Check well-formedness.\n' ex.message]) ;
            end 
            
            % store file information
            if obj.progress, waitbar(0.5,wbar,'Reading Information...'); end
            obj.XMLFileName = fileName;
            obj.FileName = binFile;
            obj.RealPath = dirName;
            obj.StudyName = studyName;
            obj.StudyPath = studyPath;
            
            obj.DOM = patdoc.getDataModelMsotProject(); % --> Whole content of the .msot XML!!
            
            obj.Name = char(obj.DOM.getScanNode.getName);
            obj.Comment = char(obj.DOM.getScanNode.getComment);
            obj.Wavelengths = zeros(length(obj.DOM.getScanNode.getWavelengths.getWavelengthArray),1,'double');
            for i = 1:length(obj.DOM.getScanNode.getWavelengths.getWavelengthArray)
                obj.Wavelengths(i) = double(obj.DOM.getScanNode.getWavelengths.getWavelengthArray(i-1).getIntValue);
            end
            
            obj.FriendlyName = char(obj.DOM.getFriendlyName);
            obj.FolderName = char(obj.DOM.getFolderName);
            obj.Scientist = char(obj.DOM.getScientist);
            obj.CreationTimeTxt = char(obj.DOM.getCreationTime.toString);
            obj.CreationTime = datenum(obj.CreationTimeTxt,'yyyy-mm-ddTHH:MM:SS.FFF');           

            obj.is3D = strcmp(char(obj.DOM.getHARDWAREDESC.getTRANSDUCER),'msot3');
            obj.is2D = strcmp(char(obj.DOM.getHARDWAREDESC.getTRANSDUCER),'msot2');
            if (obj.is3D)
                obj.TransducerType = '3D';
            else
                obj.TransducerType = '2D';
            end

            
            % presence of US data
            obj.USpresent = false;
            try obj.USpresent = logical(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDPRESENT); catch, end
            
            
            
            % Imaging Session / View FileName
            xfn = strrep(obj.FileName,'.bin','.img');
            if exist(xfn,'file')
                obj.ImgFileName = xfn;
            end
                
            if obj.progress, close(wbar); end
        end
        
        
        
        
        
        function tdom = loadImpulseResponse(obj)
            irpath = strrep((obj.XMLFileName),'msot','irf');
            if exist(irpath,'file')
                FID = fopen(irpath);
                tdom = fread(FID,'double')';
                fclose(FID);
                if obj.verbose, fprintf('Electrical Impulse Response loaded (%s)\n',irpath); end
            else 
                par.impresp = [];
                warning('Cannot load impulse response (%s), file does not exist',irpath); 
            end
        end
        
        
        % *****************************************************************
       %% Operators
       
       
       % == should compare acquisition dates (likely unique)
       function bool = eq(obj1,obj2)
           bool = (obj1.CreationTime == obj2.CreationTime);
       end
        
        
        
        
        % *****************************************************************
       %% Measurement Description
       function MD = get.MeasurementDesc(obj)
          if isempty(fieldnames(obj.MeasurementDesc))
              md = obj.DOM.getMEASUREMENTDESC;

              obj.MeasurementDesc.Averages = uint16(md.getAVERAGESPERPROJECTION); % number of shots averaged into one frame
              obj.MeasurementDesc.NumShots = uint16(md.getNUMBEROFFRAMES);       % number of shots per wavelength
              obj.MeasurementDesc.NumFrames = obj.MeasurementDesc.NumShots / obj.MeasurementDesc.Averages; % Number of frames stored in raw data
              obj.MeasurementDesc.RecordLength = uint16(md.getRECORDEDLENGTH);
              obj.MeasurementDesc.RepRate = single(md.getREPETITIONRATE.getDoubleValue);


              % wavelengths
              obj.MeasurementDesc.Wavelengths = zeros(length(md.getWAVELENGTHS.getWAVELENGTHArray),1,'double');
              for i = 1:length(md.getWAVELENGTHS.getWAVELENGTHArray)
                  obj.MeasurementDesc.Wavelengths(i) = md.getWAVELENGTHS.getWAVELENGTHArray(i-1).getIntValue;
              end
              obj.MeasurementDesc.Sequence = char(md.getSEQUENCE);

              % *** Water Absorption Coefficients
              try nwa = length(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray);
              catch, nwa = 0; end
              obj.MeasurementDesc.WaterAbsorptionCoeff = zeros(nwa,1,'single');
              for i = 1:nwa
                  obj.MeasurementDesc.WaterAbsorptionCoeff(i) = double(md.getWATERABSORPTIONCOEFFICIENTS.getWATERABSORPTIONCOEFFICIENTArray(i-1).getCoefficient);
              end
              
              % Path Length in Water (for Water correction)
              obj.MeasurementDesc.PathLengthInWater = nan;
              try obj.MeasurementDesc.PathLengthInWater = md.getPATHLENGTHINWATER.getDoubleValue; catch, end

              % *** Average Laser Energy Table (used if Average Energy Correction applied)
              try nfle = length(md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray);
              catch, nfle = 0; end
              obj.MeasurementDesc.FactoryLaserEnergyTable = zeros(nfle,1,'single');
              for i = 1:nfle
                  obj.MeasurementDesc.FactoryLaserEnergyTable(i) = single(md.getAVERAGEENERGYTABLE.getAVRAGEENERGYArray(i-1).getDoubleValue);
              end

              % Type of Laser Energy Correction
              obj.MeasurementDesc.EnergyNormalisation = char(md.getLEnergyNormalization);

              sn = obj.DOM.getScanNode;
              obj.MeasurementDesc.SWVersion = '';
              obj.MeasurementDesc.TrimSOS = 0;
              obj.MeasurementDesc.CouplantCorr = false;
              try
                  obj.MeasurementDesc.SWVersion = char(sn.getSWVersion);
                  obj.MeasurementDesc.TrimSOS = double(sn.getTrimSpeedOfSound);
                  obj.MeasurementDesc.CouplantCorr = logical(sn.getCouplantCorrection);
              catch
              end
                   
          end
          
          MD = obj.MeasurementDesc;
       end      % MeasurementDesc
       
      
       
       
       
       
       
       
       
       
       
       
       
       % ******************************************************************
      %% Hardware Description
       function HW = get.HWDesc(obj)
          if isempty(fieldnames(obj.HWDesc))
                hw = obj.DOM.getHARDWAREDESC;

                obj.HWDesc.SamplingFrequency =  hw.getSAMPLINGFREQUENCY.getDoubleValue*1e6;
                obj.HWDesc.TransducerType = char(hw.getTRANSDUCER);

                % parametric description
                obj.HWDesc.StartAngle = [];
                obj.HWDesc.StepAngle = [];
                obj.HWDesc.EndAngle = [];
                obj.HWDesc.Radius = [];
                obj.HWDesc.RadiusZ = [];
                eq = hw.getFRAMEDESC.getEQUAL;
                if (~isempty(eq))
                    obj.HWDesc.StartAngle = double(eq.getSTART);
                    obj.HWDesc.StepAngle = double(eq.getSTEP);
                    obj.HWDesc.EndAngle = double(eq.getEND);
                    obj.HWDesc.NumDetectors = round((obj.HWDesc.EndAngle - obj.HWDesc.StartAngle) / obj.HWDesc.StepAngle)+1;
                    for j = 0:numel(eq.getCONSTANTArray)-1
                       cst = eq.getCONSTANTArray(j);
                       if cst.getAxisRef == 2
                           obj.HWDesc.Radius = cst.getDoubleValue;
                       elseif cst.getAxisRef == 3
                           obj.HWDesc.RadiusZ = cst.getDoubleValue;
                       end
                    end
                end

                % sensor coordinate file
                obj.HWDesc.ProjectionData = [ ];
                try
                    pa = hw.getFRAMEDESC.getPROJECTIONArray;
                    if (~isempty(pa))
                       obj.HWDesc.NumDetectors = numel(pa);
                       obj.HWDesc.ProjectionData = zeros(numel(pa),3);
                       for j = 1:numel(pa)
                           va = pa(j).getVALUEArray;
                           % 3 coordinates
                           obj.HWDesc.ProjectionData(j,:) = [va(1).doubleValue va(2).doubleValue va(3).doubleValue ];
                           obj.HWDesc.Radius = mean(sqrt(sum(obj.HWDesc.ProjectionData.^2,2)));
                       end
                    end
                catch 
                end 

                obj.HWDesc.SpeedOfSoundBase = 0;
                try 
                    obj.HWDesc.SpeedOfSoundBase = double(hw.getSPEEDOFSOUNDBASE);
                end
                    
                obj.HWDesc.LightSpotSize = 0;
                obj.HWDesc.AxialOffset = 0;                
                obj.HWDesc.CouplantSOS = 0;
                obj.HWDesc.CouplantCOR = 0;
                obj.HWDesc.CouplantCurvature = 0;
                try
                    obj.HWDesc.LightSpotSize = double(hw.getLIGHTSPOTSIZE);
                    obj.HWDesc.AxialOffset = double(hw.getAXIALOFFSET);
                    obj.HWDesc.CouplantSOS = double(hw.getCOUPLANTSPEEDOFSOUND);
                    obj.HWDesc.CouplantCOR = double(hw.getCOUPLANTCENTEROFROTATION);
                    obj.HWDesc.CouplantCurvature = double(hw.getCOUPLANTCURVATURERADIUS);
                catch
                end
                
                sn = obj.DOM.getScanNode;
                obj.HWDesc.DeviceSN = '';
                obj.HWDesc.TransducerSN = '';
                obj.HWDesc.DAQAddress = '';
                try
                    obj.HWDesc.DeviceSN = char(sn.getDeviceSN);
                    obj.HWDesc.TransducerSN = char(sn.getTransducerSN);
                    obj.HWDesc.DAQAddress = char(sn.getDAQMACaddress);
                catch
                end
          end
          
          HW = obj.HWDesc;
       end      % hw
              
       
       
       
       
       
       
       
       
      % ******************************************************************
      %% OAMPreset
      % This is the state at time of acquisition, not necessarily the study
      % preset that is predefined
      function OA = get.OAMPreset(obj)
          if isempty(fieldnames(obj.OAMPreset))
                % ** Should go to OAM preset getter
                obj.OAMPreset = struct;
                obj.OAMPreset.Name = '';
                obj.OAMPreset.roi = [];
                try
                    oap = obj.DOM.getOAMPreset;
                    obj.OAMPreset.Name = char(oap.getName);
%                     obj.OAMPreset.Ident = char(oap.getPresetIdentifier);

                    rpa = obj.DOM.getOAMPreset.getReconPresets.getDataModelReconPreset.getSystemReconPresets;
                    rp = rpa(1).getReconstructionPreset;
                    obj.OAMPreset.RoiHigh = double(rp.getRoiHigh);
                    obj.OAMPreset.RoiLow = double(rp.getRoiLow);
                    obj.OAMPreset.roi = double(rp.getRoiHigh) - double(rp.getRoiLow);
                    obj.OAMPreset.TimeRes = rp.getTimeRes;
                    obj.OAMPreset.Projections = rp.getProjections;
                    obj.OAMPreset.n = rp.getResolution;
                    
                    rpm = obj.DOM.getOAMPreset.getReconPresets.getDataModelReconPreset;
                    obj.OAMPreset.FilterLow = rpm.getFilterLow;
                    obj.OAMPreset.FilterHigh = rpm.getFilterHigh;
                    obj.OAMPreset.DepthCorrection = logical(rpm.getDepthCorrection);
                    obj.OAMPreset.BackgroundAbsorption = double(rpm.getBackgroundAbsorption);
                    obj.OAMPreset.BackgroundOxygenation = rpm.getBackgroundOxygenation;
                    obj.OAMPreset.UserSoundTrim = rpm.getUserSoundTrim;
                    
                    mp = obj.DOM.getOAMPreset.getMspPresets.getDataModelMspPreset;
                    obj.OAMPreset.MSP.DiscardNegatives = mp.getDiscardNegativeValues;
                    obj.OAMPreset.MSP.Method = char(mp.getMethod);
                    obj.OAMPreset.MSP.BGWavelength = mp.getBgWavelength.getIntValue;
                    obj.OAMPreset.MSP.Spectra = cell(mp.getUserSelectedSpectra.getStringArray);
                    
                    sett = obj.DOM.getOAMPreset.getImagingSettingsPreset;
                    obj.OAMPreset.ViewSettings.UltrasoundMinimumScaling = double(sett.getUltrasoundMinimumScaling);
                    obj.OAMPreset.ViewSettings.UltrasoundMaximumScaling = double(sett.getUltrasoundMaximumScaling);
                    obj.OAMPreset.ViewSettings.BackgroundMinimumScaling = double(sett.getBackgroundMinimumScaling);
                    obj.OAMPreset.ViewSettings.BackgroundMaximumScaling = double(sett.getBackgroundMaximumScaling);
                    obj.OAMPreset.ViewSettings.ForegroundMinimumScaling = double(sett.getForegroundMinimumScaling);
                    obj.OAMPreset.ViewSettings.ForegroundMaximumScaling = double(sett.getForegroundMaximumScaling);
                    obj.OAMPreset.ViewSettings.Visible3DGridPlanesTypes = char(sett.getVisible3DGridPlanesTypes);

                catch
                end
          end
          
          OA = obj.OAMPreset;
      end
      
       
       
       
       
       % ******************************************************************
        %% US Information
        function US = get.US(obj)
            if ~obj.USpresent, US = []; return; end
            if ~obj.structureLoaded, calculateStructure(obj); end
            
            if isempty(fieldnames(obj.US))
                % US Info (might not be present)
                obj.US.timestamps = single([]);
                obj.US.OAframes = uint32([]);
                obj.US.OArep = uint16([]);

                % resolution
                obj.US.n = 0;
                obj.US.nz = 0;
                try obj.US.n = uint16(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDRESOLUTION); catch, end
                try obj.US.n = uint16(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDFIELDOFVIEW.getPixelCount.getX); catch, end
                try obj.US.nz = uint16(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDFIELDOFVIEW.getPixelCount.getZ); catch, end
                
                % field of view
                obj.US.roi = 0.0;
                obj.US.roiz = 0.0;
                try obj.US.roi = double(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDFIELDOFVIEW.getExtents.getX); catch, end
                try obj.US.roiz = double(obj.DOM.getMEASUREMENTDESC.getULTRASOUNDFIELDOFVIEW.getExtents.getZ); catch, end
                                                
                for i = 1:numel(obj.ScanFrames)
                    uoff = obj.ScanFrames(i).USOffset;
                    if numel(obj.US.timestamps) < uoff+1 
                        obj.US.timestamps(uoff+1) = obj.ScanFrames(i).RelTime;
                        obj.US.OAframes(uoff+1) = uint32(i);
                        obj.US.OArep(uoff+1) = uint16(obj.ScanFrames(i).Repetition);
                    end
                end
            end
            
            US = obj.US;
        end
       
       
       
       
       
        
        
       
       % ******************************************************************
      %% Structural Information (requires complete parse)
       function SF = get.ScanFrames(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SF = obj.ScanFrames;
       end
       
       function SF = get.ZPositions(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SF = obj.ZPositions;
       end
      
       function RN = get.RepNum(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          RN = obj.RepNum;
       end
     
       function RN = get.RunNum(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          RN = obj.RunNum;
       end

       function RN = get.ShotNum(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          RN = obj.ShotNum;
       end

       function AT = get.AverageTemperature(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          AT = obj.AverageTemperature;
       end
       
       function SS = get.ScanStructure(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.ScanStructure;
       end
       
       function SS = get.RelTime(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.RelTime;
       end
       
       function SS = get.LaserEnergy(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.LaserEnergy;
       end
       
       function SS = get.DiodeReadout(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.DiodeReadout;
       end
       
       function SS = get.CorrectionFactor(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.CorrectionFactor;
       end
       
       function SS = get.Temperature(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = obj.Temperature;
       end
       
       function SS = get.Timestamps(obj)
          if ~obj.structureLoaded, calculateStructure(obj); end
          SS = squeeze(obj.RelTime(:,1,1,1,1));
       end
              

       
              
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       % ******************************************************************
      %% Structure Parsing
       
       function calculateStructure(obj)
            tic
            if obj.progress, wbar = waitbar(0,'Reading Structural Information...'); end

            obj.structureLoaded = true; % do this now to avoid recursion
            
            try
                sn = obj.DOM.getScanNode;
                fa = sn.getScanFrames.getDataModelScanFrameArray;

                zpos = single([]);
                ts0 = double(0);
                nfr = length(fa);
                obj.ScanFrames(1) = scanFrame();  % initialise array
                for i = 1:nfr
                    if obj.progress && (mod(i,25) == 1), wbar = waitbar(i/nfr*0.75,wbar,'Reading Frame Information...'); end
                    fr = fa(i).getFrame; 

                    obj.ScanFrames(i).IDOffset = uint32(fa(i).getIDOffset);
                    obj.ScanFrames(i).meta = obj;
                    obj.ScanFrames(i).FrameNum = uint32(i);

                    try wl = double(fr.getWAVELENGTH.getIntValue); catch, wl = double(nan); end
                    obj.ScanFrames(i).Wavelength = wl;
                    try zr = single(fr.getZPOS.getDoubleValue); catch, zr = single(0); end
                    obj.ScanFrames(i).ZPosReal = zr; 

                    % Check if the POWER field exists and contains a value. If not, try to retrieve value from the OverallCorrectionFactor field
                    try 
                        obj.ScanFrames(i).CorrectionFactor = single(fr.getPOWER.getDoubleValue); 
                    catch
                        try 
                            obj.ScanFrames(i).CorrectionFactor = single(double(fr.getOverallCorrectionFactor)); 
                            % OverallCorrectionFactor field is not a 'complexType', thus fetching the value is performed without supplementary method getDoubleValue 
                        catch
                            obj.ScanFrames(i).CorrectionFactor = single(nan); 
                        end
                    end
                    
                    try obj.ScanFrames(i).Temperature = single(fr.getTEMPERATURE.getDoubleValue); catch,  obj.ScanFrames(i).Temperature = single(0); end  

                    ts = fr.getTimestamp.doubleValue * 1e-7;
                    if (i == 1), ts0 = ts; end
                    obj.ScanFrames(i).RelTime = double(ts - ts0);
                    obj.ScanFrames(i).Timestamp = double(ts);

                    % Later Additions (might not be present)
                    obj.ScanFrames(i).DiodeReadout = single(nan);
                    obj.ScanFrames(i).LaserEnergy = single(nan);
                    try obj.ScanFrames(i).LaserEnergy = single(fr.getLASERENERGY.getDoubleValue); catch, end
                    try obj.ScanFrames(i).DiodeReadout = single(fr.getDIODEREADOUT.getDoubleValue); catch, end

                    obj.ScanFrames(i).USOffset = uint32(nan);
                    try
                        uoff = uint32(fr.getUltraSoundFrameOffset);
                        obj.ScanFrames(i).USOffset = uoff;
                    catch
                    end

                    ru = uint16(fr.getRUN);
                    if ru == 0, ru = uint16(1); end           % if no RUNs, then it will be zero
                    obj.ScanFrames(i).Run = ru;
                    obj.ScanFrames(i).Repetition = uint16(fr.getREPETITION);

                    % Find ZPositions (with Rounding Problems)
                    [m, zp] = min(abs(zpos - zr));
                    if isempty(m) || m > 0.05 
                        zp = numel(zpos)+1;
                        zpos(zp) = zr;
                    end
                    obj.ScanFrames(i).ZNum = uint16(zp);

                    % Wavelength
                    obj.ScanFrames(i).WLNum = uint16(find(obj.Wavelengths == wl));
                end
                obj.RepNum = max([obj.ScanFrames.Repetition]);
                obj.RunNum = max([obj.ScanFrames.Run]); 
                obj.AverageTemperature = mean([obj.ScanFrames.Temperature]);

                % round and assign ZPositions similar to ViewMSOT
                zr = [obj.ScanFrames.ZPosReal];
                zn = [obj.ScanFrames.ZNum];
                for jj = 1:numel(zpos)
                   zp = round(mean(zr(zn == jj))*10)/10;
                   zr(zn == jj) = zp;
                   [obj.ScanFrames(zn == jj).ZPos] = deal(zp);
                   zpos(jj) = zp;
                end
                obj.ZPositions = zpos;

                % Number of Shots
                numshots = obj.MeasurementDesc.NumFrames;
                if (numshots == 0) || (numshots > 1 && numshots*(obj.RunNum*numel(obj.ZPositions)*obj.RepNum*numel(obj.Wavelengths)) > numel(obj.ScanFrames)) 
                   % Legacy support for datasets without Number of Frames
                   % support for datasets with NumFrames not set
                   % accordingy (like vM 4.0.0.12)
                   numshots = round(numel(obj.ScanFrames)./(obj.RunNum*numel(obj.ZPositions)*obj.RepNum*numel(obj.Wavelengths)));
                   obj.MeasurementDesc.NumFrames = numshots;
                end
                obj.ShotNum = numshots;

                % Initialise Structure Matrices
                obj.ScanStructure = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),numshots,'single');
                obj.RelTime = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'double');
                obj.LaserEnergy = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                obj.DiodeReadout = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                obj.CorrectionFactor = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                obj.Temperature = nan(obj.RunNum,numel(obj.ZPositions),obj.RepNum,numel(obj.Wavelengths),'single');
                obj.idoffset = zeros(max([obj.ScanFrames.IDOffset])+1,1,'uint32');

                oru = 0; ore = 0; ozp = 0; owl = 0; sh = 1;
                for i = 1:nfr
                    if obj.progress && (mod(i,25) == 1), wbar = waitbar(0.75+i/nfr*0.25,wbar,'Reading Structural Information...'); end
                    fr = obj.ScanFrames(i);

                    % reverse linking for recons and msps
                    obj.idoffset(obj.ScanFrames(i).IDOffset+1) = i;  % map id offset to frame number

                    % shotnumber
                    if (oru ~= fr.Run) || (ore ~= fr.Repetition) || (ozp ~= fr.ZNum) || (owl ~= fr.WLNum)
                        sh = 1;
                        obj.RelTime(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum) = fr.RelTime;
                        obj.Temperature(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum) = fr.Temperature;
                    else
                        sh = sh + 1;
                    end
                    obj.ScanFrames(i).ShotNum = uint16(sh);

                    obj.ScanStructure(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum,obj.ScanFrames(i).ShotNum) = i;
                    obj.LaserEnergy(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum,obj.ScanFrames(i).ShotNum) = fr.LaserEnergy;
                    obj.DiodeReadout(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum,obj.ScanFrames(i).ShotNum) = fr.DiodeReadout;
                    obj.CorrectionFactor(fr.Run,fr.ZNum,fr.Repetition,fr.WLNum,obj.ScanFrames(i).ShotNum) = fr.CorrectionFactor;

                    % memorise
                    oru = fr.Run; ore = fr.Repetition; ozp = fr.ZNum; owl = fr.WLNum; 
                end

                obj.structureTime = single(toc);
                if obj.verbose, fprintf('Structure Information read in %.1fs\n',obj.structureTime); end
                if obj.progress, close(wbar); end
                
               
            catch ex
                if obj.progress, close(wbar); end
                obj.structureLoaded = false;
                rethrow(ex);
            end
       end
 
       
       
       
 
       
       
       
       
       
       
       
       
       
       
       
       
       
       % ******************************************************************
      %% View Information      
       function IS = get.ImgSession(obj)
          if isempty(obj.ImgSession) && ~isempty(obj.ImgFileName)
              try
                  obj.ImgSession = loadMSOTView(obj.ImgFileName,obj);
              end
          end
          
          IS = obj.ImgSession;
       end
      
       
       
       
       

       
       
       
       
       
       
       
       
       % ******************************************************************
      %% Recon Node Information      
       function RN = get.ReconNode(obj)
          if isempty(obj.ReconNode)
                try 
                    rna = obj.DOM.getReconNodes.getDataModelNewReconstructionNodeArray;
                catch
                    RN = []; 
                    return 
                end

                if obj.progress 
                    wbar = waitbar(0,'Reading Recon Information...'); 
                end

                % read all Recons
                for rn = 1:length(rna)
                    if obj.progress 
                        wbar = waitbar(rn/length(rna),wbar,'Reading Recon Information...'); 
                    end
                    
                    obj.ReconNode(rn) = reconNode(obj,rna(rn),uint16(rn));
                end
        
                if obj.progress 
                    close(wbar); 
                end
          end
          
          RN = obj.ReconNode;
       end
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       % ******************************************************************
      %% MSP Node Information      
       function MN = get.MSPNode(obj)
          if isempty(obj.MSPNode)
                try mna = obj.DOM.getMspNodes.getDataModelNewMspNodeArray;
                catch, MN = []; return; end

                if obj.progress, wbar = waitbar(0,'Reading MSP Information...'); end

                for mn = 1:length(mna)
                    if obj.progress, waitbar(mn/length(mna),wbar) ; end
                    obj.MSPNode(mn) = mspNode(obj,mna(mn),uint16(mn));
                end
                
                if obj.progress, close(wbar); end
          end
          
          MN = obj.MSPNode;
       end
       
       
       
       
       
       
       
       % ******************************************************************
      %% String Representations
      function c = get.RunString(obj)
          run = (1:obj.RunNum)';
          ts = squeeze(obj.RelTime(:,1,1,1,1));
          hour = floor(ts/60/60);
          min = floor(mod(ts,60*60)/60);
          sec = floor(mod(ts,60));

          if numel(run) ~= numel(ts), warning('Unsupported Dataset, number of timestamps != number of repetitions'); c = {}; end
          str = [num2str(hour,'%02i:') num2str(min,'%02i:') num2str(sec,'%02i') num2str(run,'x(%03i)')];
          c = mat2cell(str,ones(numel(ts),1));
          c = regexprep(c,'x',' ');
      end
      
      function c = get.ZPositionString(obj)
          if isempty(obj.ZPositions), c = {}; return; end 
          c = mat2cell(num2str(reshape(obj.ZPositions,numel(obj.ZPositions),1),'%.1fmm'),ones(numel(obj.ZPositions),1));
      end
      
      function c = get.RepString(obj)
          reps = (1:obj.RepNum)';
          ts = squeeze(obj.RelTime(1,1,:,1,1));
          min = floor(ts/60);
          sec = floor(mod(ts,60));
          msec = floor(mod(ts,1)*100);

          if numel(reps) ~= numel(ts), warning('Unsupported Dataset, number of timestamps != number of repetitions'); c = {}; end
          str = [num2str(min,'%02i:') num2str(sec,'%02i.')  num2str(msec,'%02i') num2str(reps,'x(%03i)')];
          c = mat2cell(str,ones(numel(ts),1));
          c = regexprep(c,'x',' ');
      end
      
      function c = get.WavelengthString(obj)
          if isempty(obj.Wavelengths), c = {}; return; end 
          c = mat2cell(num2str(obj.Wavelengths,'%inm'),ones(numel(obj.Wavelengths),1));
      end
       
       

      
      
       % ******************************************************************
      %% Laser Energy Plot
      function plotEnergy(obj,plotDiode)
        if nargin == 1, plotDiode = false; end
        
        yt = 10;
        yt2 = 500;
        ntick = ceil(max(obj.LaserEnergy(:))/yt);
        
        figure;
        if plotDiode
            [ax, h1, h2] = plotyy(obj.Wavelengths, squeeze(mean(obj.LaserEnergy,2)), obj.Wavelengths, squeeze(mean(obj.DiodeReadout,2)));
            h2.Marker = '.';
            xlim(ax(2),[min(obj.Wavelengths) max(obj.Wavelengths)]);
            ylim(ax(2),[0 ntick*yt2]);
            ylabel(ax(2),'Diode Readout (a.u.)');
            set(ax(2),'YTick',0:yt2:ntick*yt2);
        else
            h1 = plot(obj.Wavelengths, squeeze(mean(obj.LaserEnergy,2)));
            ax = gca;
        end
        
        xlabel(ax(1),'Wavelengths (nm)');
        ylabel(ax(1),'Laser Energy (mJ)');
        xlim(ax(1),[min(obj.Wavelengths) max(obj.Wavelengths)]);
        ylim(ax(1),[0 ntick*yt]);
        h1.Marker = '.';
        set(ax(1),'YTick',0:yt:ntick*yt);
        set(ax(1),'XTick',0:50:4000);
        title(['Integrated Laser Energy Meter (' obj.Name ')']);
        grid on;      
      end

    end       % methods
end      % classdef