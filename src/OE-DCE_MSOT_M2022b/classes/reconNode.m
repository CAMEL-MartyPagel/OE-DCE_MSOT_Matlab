classdef reconNode < handle
    properties
        meta msotData;
        nodeNumber uint16;
        DataModelVersion char;
        
        Name char;
        Comment char;
        GUID char;
        Method char;
        Resolution uint16;
        ResolutionY uint16;
        ZResolution uint16;
        Projections uint16;
        ROI double;

        % Dependent on structure calculation
        Frames struct;
        ZPositions single;
        ZNum uint16;
        Runs uint16;
        Repetitions uint16;
        Wavelengths double;
        
        ReconStructure uint32;
        RelTime single ;
        Timestamps single ;
        IDLookup uint32 ;
        
        FRAMEARRAY;
        DOM;
    end
    
    properties ( Access = private )
        structureLoaded logical = false;
    end
    
    methods
        
       % ******************************************************************
      %% Constructor
       function obj = reconNode(varargin)
           if nargin ~= 3, return; end
           
           obj.meta = varargin{1};
           obj.DOM = varargin{2};
           obj.nodeNumber = varargin{3};
           obj.DataModelVersion = char(obj.DOM.getDataModelVersion);
           
           obj.Name = char(obj.DOM.getName);
           obj.Comment = char(obj.DOM.getComment);
           obj.GUID = char(obj.DOM.getGUID);
           obj.Method = char(obj.DOM.getMethod);
           obj.Resolution = uint16(obj.DOM.getResolution);
           obj.ResolutionY = obj.Resolution;
           obj.ZResolution = uint16(obj.DOM.getNumOfFrames);
           obj.ROI = double(obj.DOM.getRoi);

           if obj.Resolution == 0 % since FoV was introduced
               fov = obj.DOM.getFIELDOFVIEW.getPixelCount;
               obj.Resolution = uint16(fov.getX);
               obj.ResolutionY = uint16(fov.getY);
               obj.ZResolution = uint16(fov.getZ);
               % for 2D now the orientation in ViewMSOT is (x,z)
               if obj.ResolutionY == 1
                   obj.ResolutionY = obj.ZResolution;
                   obj.ZResolution = uint16(1);
               end
           end
           if isempty(obj.ROI)
               ext = obj.DOM.getFIELDOFVIEW.getExtents;
               obj.ROI = double(ext.getX);
           end
           
           obj.Projections = uint16(obj.DOM.getProjections); 
           obj.FRAMEARRAY = obj.DOM.getReconFrames.getDataModelReconstructionFrameArray;           
       end
        
       
       
       
       
       
       % ******************************************************************
      %% Dependent Accessors
        function FA = get.Frames(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Frames;
        end
        
        function FA = get.Runs(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Runs;
        end
        function FA = get.Repetitions(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Repetitions;
        end
        function FA = get.ZNum(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.ZNum;
        end
        function FA = get.Wavelengths(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Wavelengths;
        end
        function FA = get.ZPositions(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.ZPositions;
        end
        function FA = get.ReconStructure(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.ReconStructure;
        end
        function FA = get.RelTime(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.RelTime;
        end
        function FA = get.Timestamps(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Timestamps;
        end
        function ID = get.IDLookup(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            ID = obj.IDLookup;
        end

        
        
        
        
        
        
        
        
       % ******************************************************************
      %% Structure Parsing
        function obj = calculateStructure(obj)
            if obj.meta.progress, wbar = waitbar(0,'Reading Recon Frames...'); end
            obj.structureLoaded = true;      % avoid recursion

            try
                sv = size(obj.meta.ScanStructure);
                if numel(sv) < 4, sv = [sv ones(1,4-numel(sv))]; end
                nfra = length(obj.FRAMEARRAY);

                idlookup = zeros(nfra*obj.ZResolution,1,'uint32');
                ridm = zeros(size(obj.meta.ScanStructure),'uint32');
                tsm = zeros(size(obj.meta.ScanStructure),'single');
                averaged = false();
                for rfr = 1:nfra
                    fr = obj.FRAMEARRAY(rfr);
                    if obj.meta.progress && (mod(rfr,25) == 1), wbar = waitbar(rfr/nfra,wbar,'Reading Recon Frames...'); end
                    obj.Frames(rfr).IDOffset = fr.getIDOffset;
                    obj.Frames(rfr).ScanFrames = str2num(fr.getScanFramesID.getStringValue);     
                    % use wavelength and zpos from first frame
                    frid = obj.meta.idoffset(obj.Frames(rfr).ScanFrames(1)+1);
                    obj.Frames(rfr).Wavelength = obj.meta.ScanFrames(frid).Wavelength;
                    obj.Frames(rfr).ZPos = obj.meta.ScanFrames(frid).ZPos;
                    obj.Frames(rfr).ZNum = obj.meta.ScanFrames(frid).ZNum;
                    obj.Frames(rfr).ShotNum = obj.meta.ScanFrames(frid).ShotNum;
                    obj.Frames(rfr).Run = obj.meta.ScanFrames(frid).Run;
                    obj.Frames(rfr).Repetition = obj.meta.ScanFrames(frid).Repetition;
                    obj.Frames(rfr).Timestamp = obj.meta.ScanFrames(frid).Timestamp;
                    obj.Frames(rfr).RelTime = obj.meta.ScanFrames(frid).RelTime;

                    % get all power values
                    for j = 1:length(obj.Frames(rfr).ScanFrames)
                        obj.Frames(rfr).LaserEnergy(j) = obj.meta.ScanFrames(obj.meta.idoffset(obj.Frames(rfr).ScanFrames(j)+1)).LaserEnergy;
                    end
                    if ~averaged && length(obj.Frames(rfr).LaserEnergy) > 1
                        averaged = true();
                    end

                    % assign id to correct entry in structure
                    ru = find(obj.Runs == obj.Frames(rfr).Run);
                    if isempty(ru) 
                        ru = numel(obj.Runs)+1; 
                        obj.Runs = [obj.Runs, obj.Frames(rfr).Run]; 
                    end
                    zp = find(obj.ZNum == obj.Frames(rfr).ZNum);
                    if isempty(zp)
                        zp = numel(obj.ZNum)+1; 
                        obj.ZNum = [obj.ZNum obj.Frames(rfr).ZNum]; 
                    end
                    re = find(obj.Repetitions == obj.Frames(rfr).Repetition);
                    if isempty(re)
                        re = numel(obj.Repetitions)+1; 
                        obj.Repetitions = [obj.Repetitions obj.Frames(rfr).Repetition]; 
                    end
                    wl = find(obj.Wavelengths == obj.Frames(rfr).Wavelength);
                    if isempty(wl) 
                        wl = numel(obj.Wavelengths)+1; 
                        obj.Wavelengths = [obj.Wavelengths obj.Frames(rfr).Wavelength]; 
                    end
                    sh = obj.meta.ScanFrames(frid).ShotNum;
                    ridm(ru,zp,re,wl,sh) = rfr;
                    tsm(ru,zp,re,wl,sh) = obj.Frames(rfr).RelTime;
                    idlookup(obj.Frames(rfr).IDOffset+1) = rfr;    % reverse lookup to find when loading MSP in 3D cases, where IDOffset != Frame Number
                    clear ru zp re wl sh;
                end
                
                if averaged  % remove average dimension
                    ridm = ridm(:,:,:,:,1);
                end
                
                % save ReconStructure as matrix
                obj.ZPositions = obj.meta.ZPositions(obj.ZNum);
                obj.ReconStructure = ridm;
                obj.RelTime = tsm;
                obj.Timestamps = squeeze(obj.RelTime(:,1,1,1,1));
                obj.IDLookup = idlookup;
                clear ridm tsm ridc;
                if obj.meta.progress, close(wbar); end
                
            catch ex
                if obj.meta.progress, close(wbar); end
                obj.structureLoaded = false;
                rethrow(ex);
            end

        end
        
    end
end