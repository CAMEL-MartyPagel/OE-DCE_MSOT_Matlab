classdef mspNode < handle,
    properties
        meta@msotData;
        nodeNumber@uint16;
        DataModelVersion@char vector;
        
        Name@char vector;
        Comment@char vector;
        GUID@char vector;
        Method@char vector;
        InputSpectra@cell vector;

        ReconNode@reconNode;
        ReconGUID@char vector;
        ReconNodeID@uint16;
        ReconMethod@char vector;
        Resolution@uint16;
        ResolutionY@uint16;
        ZResolution@uint16;
        Projections@uint16;
        ROI@double;

        % Dependent on structure calculation
        Slices@struct vector;
        ZPositions@single vector;
        ZNum@uint16 vector;
        Runs@uint16 vector;
        Repetitions@uint16 vector;
        Wavelengths@double vector;
        
        Structure@uint32;
        RelTime@single ;
        Timestamps@single ;
        SliceIndex@uint32 ;
        
        FRAMEARRAY;
        DOM;       
    end
    
    properties ( Access = private )
        structureLoaded@logical = false;
    end
    
    methods
        
       % ******************************************************************
      %% Constructor
       function obj = mspNode(varargin)
           if nargin ~= 3, return; end
           
           obj.meta = varargin{1};
           obj.DOM = varargin{2};
           obj.nodeNumber = varargin{3};
           
            obj.DataModelVersion = char(obj.DOM.getDataModelVersion);
            obj.Name = char(obj.DOM.getName);
            obj.Comment = char(obj.DOM.getComment);
            obj.Method = char(obj.DOM.getMethod);
            obj.GUID = char(obj.DOM.getGUID);
            
            obj.ReconGUID = char(obj.DOM.getReconGUID);
            % Related Wavelengths
            wla = obj.DOM.getRelatedWavelengths.getWavelengthArray;
            for wl = 1:numel(wla)
                obj.Wavelengths(wl) = double(wla(wl).getIntValue);
            end
            % Input Spectra
            obj.InputSpectra = cell(obj.DOM.getInputSpectra.getStringArray);
            
            % Find Recon Info
            rn = 1;
            while ~strcmp(obj.meta.ReconNode(rn).GUID,obj.ReconGUID)
                rn = rn + 1;
                if (rn > numel(obj.meta.ReconNode)), rn = 0; break; end;
            end
            obj.ReconNodeID = uint16(rn);
            if (rn > 0)
                obj.ReconNode = obj.meta.ReconNode(rn);
                obj.ReconMethod = obj.meta.ReconNode(rn).Method;
                obj.Resolution = obj.meta.ReconNode(rn).Resolution;
                obj.ResolutionY = obj.meta.ReconNode(rn).ResolutionY;
                obj.ZResolution = obj.meta.ReconNode(rn).ZResolution;
                obj.Projections = obj.meta.ReconNode(rn).Projections;
                obj.ROI = obj.meta.ReconNode(rn).ROI;
            else
                warning(['Recon ' obj.ReconGUID ' not found.\n']);
            end

            obj.FRAMEARRAY = obj.DOM.getSlicesFrame.getMspSliceFrameArray;
       end
    
    
    
    
    
    
    
       % ******************************************************************
      %% Dependent Accessors
        function FA = get.Slices(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Slices;
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

        function FA = get.ZPositions(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.ZPositions;
        end
        function FA = get.Structure(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Structure;
        end
        function FA = get.RelTime(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.RelTime;
        end
        function FA = get.Timestamps(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            FA = obj.Timestamps;
        end
        function ID = get.SliceIndex(obj)
            if ~obj.structureLoaded, calculateStructure(obj); end
            ID = obj.SliceIndex;
        end

        
        
        
        
        
        
        
       % ******************************************************************
      %% Structure Parsing
        function obj = calculateStructure(obj)
            if obj.meta.progress, wbar = waitbar(0,'Reading MSP Frames...'); end;
            obj.structureLoaded = true;      % avoid recursion

            try
                nmfr = numel(obj.FRAMEARRAY);
                for mfr = 1:nmfr
                    if obj.meta.progress && (mod(mfr,25) == 1), waitbar(mfr/nmfr,wbar,'Reading MSP Frames...'); end;
                    fr = obj.FRAMEARRAY(mfr);

                    obj.Slices(mfr).HasErrors = logical(fr.getHasErrors);
                    % Recon Info
                    if obj.ReconNodeID > 0,
                        obj.Slices(mfr).ReconFrames = str2num(fr.getReconFrameID.getStringValue);
                        % reverse lookup to find recon frame in 3D cases, where IDOffset != Frame Number
                        rframes = obj.ReconNode.IDLookup(obj.Slices(mfr).ReconFrames(:)+1);  % all associated recon frames
                        sframeid = obj.ReconNode.Frames(rframes(1)).ScanFrames(1)+1; % first scan frame for referencing
                        obj.Slices(mfr).ReconFrameNumbers = rframes;
                        obj.Slices(mfr).ScanFrames = cell2mat({obj.ReconNode.Frames(rframes).ScanFrames});
                        obj.Slices(mfr).ZPos = obj.meta.ScanFrames(sframeid).ZPos;
                        obj.Slices(mfr).ZNum = obj.meta.ScanFrames(sframeid).ZNum;
                        obj.Slices(mfr).Run = obj.meta.ScanFrames(sframeid).Run;
                        obj.Slices(mfr).Repetition = obj.meta.ScanFrames(sframeid).Repetition;
                        obj.Slices(mfr).RelTime = obj.meta.ScanFrames(sframeid).RelTime;
                        zp = find(obj.ZNum ==  obj.Slices(mfr).ZNum);
                        % if z-position was not yet added, add it now and keep the index of
                        % the added z-positions for the structure indexing
                        if (isempty(zp)) 
                            obj.ZNum = [obj.ZNum obj.Slices(mfr).ZNum];
                            zp = numel(obj.ZNum);
                        end
                    end
                    % Component Info
                    ca = fr.getComponentList.getMspSliceComponentArray;
                    for c = 1:numel(ca)
                       obj.Slices(mfr).Components(c).IDOffset = ca(c).getIDOffset;
                       obj.Slices(mfr).Components(c).ComponName = char(ca(c).getComponName);
                       obj.Slices(mfr).Components(c).ComponRef = char(ca(c).getComponRef);
                       obj.Slices(mfr).Components(c).Invert = logical(ca(c).getInvert);
                       obj.Slices(mfr).Components(c).SpectrumName = char(ca(c).getSpectraName);
                       mstr(obj.Slices(mfr).Run,zp,obj.Slices(mfr).Repetition,c) = ca(c).getIDOffset+1;
                       slref(ca(c).getIDOffset+1) = mfr;
                    end
                end
                obj.ZPositions = obj.meta.ZPositions(obj.ZNum);
                obj.Runs = obj.ReconNode.Runs;
                obj.Repetitions = obj.ReconNode.Repetitions;
                obj.Timestamps = obj.ReconNode.Timestamps;
                obj.RelTime = obj.ReconNode.RelTime;
                obj.Structure = uint32(mstr);
                obj.SliceIndex = uint32(slref);                
                
                if obj.meta.progress, close(wbar); end;
            catch ex,
                if obj.meta.progress, close(wbar); end;
                obj.structureLoaded = false;
                rethrow(ex);
            end

        end        
        
        
        
        
        
    end  
    
end