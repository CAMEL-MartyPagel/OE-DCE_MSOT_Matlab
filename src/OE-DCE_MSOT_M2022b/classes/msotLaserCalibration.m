classdef msotLaserCalibration < handle
    properties
        FileName@char vector = char([]);
        FileVersion@uint16 scalar = uint16(1);
        DescriptiveHeaderText@char vector = char([]);
        SafetyNumber@uint16 scalar = uint16(0);

        PulseDividerCount@uint16 scalar = uint16(1);
        CoefficientCount@uint16 scalar = uint16(6);
        
        LaserUndividedTrigger@uint16 scalar = uint16(100);
        PulseDividerValues@uint16 vector = uint16([]);
        Wavelengths@uint16 vector = uint16([]);
        
        Elements@msotLaserCalibrationElement vector = msotLaserCalibrationElement.empty(0,1);   % Calibration Elements (one per WL and Divider)
    end
    
    properties (SetAccess = private)
     
    end    
    
    methods 
        function obj = msotLaserCalibration(FileName)
            % Constructor (takes FileName as optional argument)
            if nargin == 0, return; end
            
            if nargin >= 1,
                readFile(obj,FileName);
            end
        end
        
        function Value = Bytes2Uint16(~,FirstSecondByte)
            Value = FirstSecondByte(1) + FirstSecondByte(2) * 256;
        end

        function readFile(obj,FileName)
            % Populates Class Values from Calibration File
            [fid, errmsg] = fopen(FileName ,'r');
            if fid < 0,
                error('Error opening Laser Calibration File: %s',errmsg);
            end

            % Reading at once
            DataBytes=fread(fid);
            fclose(fid);

            obj.FileName = FileName;
            obj.FileVersion = obj.Bytes2Uint16(DataBytes(1:2));

            DescriptiveHeaderCount= obj.Bytes2Uint16(DataBytes(3:4));
            obj.DescriptiveHeaderText = native2unicode(DataBytes(5: (DescriptiveHeaderCount+4)),'UTF-8')';
            Count= (DescriptiveHeaderCount +5) ;
            
            obj.SafetyNumber = obj.Bytes2Uint16(DataBytes(Count:(Count+1)));Count = Count +2;
            obj.LaserUndividedTrigger = obj.Bytes2Uint16(DataBytes(Count:(Count+1)));Count = Count + 2 ;
            obj.PulseDividerCount = obj.Bytes2Uint16(DataBytes(Count:(Count+1)));Count = Count + 2 ;
            obj.CoefficientCount = obj.Bytes2Uint16(DataBytes(Count:(Count+1)));Count = Count + 2;

            obj.PulseDividerValues = zeros(obj.PulseDividerCount,1,'uint16');
            for i = 1: obj.PulseDividerCount
                obj.PulseDividerValues(i) = obj.Bytes2Uint16(DataBytes(Count:(Count+1)));
                Count = Count + 2;
            end

            wlc = obj.Bytes2Uint16(DataBytes(Count:(Count+1))); Count = Count + 2;
            obj.Wavelengths = zeros(wlc,1,'uint16');
            cc = 0;
            for k = 1 : wlc
                wl = obj.Bytes2Uint16(DataBytes(Count:(Count+1))); Count = Count + 2;
                obj.Wavelengths(k) = wl;
                for l = 1 : obj.PulseDividerCount
                    dp = obj.PulseDividerValues(l);
                    ae = typecast(uint8(DataBytes(Count:(Count+7))),'double'); Count = Count + 8;
                    att = typecast(uint8(DataBytes(Count:(Count+7))),'double'); Count = Count + 8;
                    co = zeros(obj.CoefficientCount ,1);
                    for  m = 1 : obj.CoefficientCount 
                        co(m) = typecast(uint8(DataBytes(Count:(Count+7))),'double'); Count = Count + 8;
                    end
                    cc = cc+1;
                    obj.Elements(cc) = msotLaserCalibrationElement(wl,dp,ae,att,co);
                end
            end
        end
        
        
        function lst = getElements(obj,wls,pd)
            % Gets a list of Elements based on the specified conditions:
            % 1: Wavelengths (list of wavelenghts, all if empty)
            % 2: Pulse Divider Values (list, all if empty)
            if nargin <= 1 || isempty(wls)
                wls = obj.Wavelengths;
            end
            if nargin <= 2 || isempty(pd),
                pd = obj.PulseDividerValues(1);
            end
            
            % build 
            loc = ismember([obj.Elements.Wavelength],wls) & ismember([obj.Elements.PulseDivider],pd);
            lst = obj.Elements(loc);
        end
        
        function en = calculate(obj,in1,in2)
            % Calculate Energy based on Diode Readout, either takes a
            % msotFrame or Wavelength and DiodeReadout as parameters
            if nargin == 2 && isa(in1,'msotFrame'),
                el = obj.getElements(in1.Wavelength);
                en = el.calculate(in1.DiodeReadout);
            else
                el = obj.getElements(in1);
                en = el.calculate(in2);
            end
        end
    end

end