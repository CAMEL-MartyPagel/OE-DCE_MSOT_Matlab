classdef (Abstract) msotProbe < msotSettings
    % msotProbe Abstract class for MSOT Probe Definition
    %   Contains all common Probe Properties, derived classes available for
    %   Parametric and Coordinate Definitions exist:
    %   - msotProbeCoordinates
    %   - msotProbeParametric
    %
    % msotProbe Properties:
    % *META Properties*
    %   ProbeID - iThera Probe ID - mandatory (e.g. PR010 - see approriate List)
    %   DesignID - Detector Design ID
    
properties
    ProbeID = 'unknownProbe'; % iThera Probe ID - mandatory (e.g. PR010 - see approriate List) --> PType
    DesignID = 'cdc12358-1';  % Detector Design ID
end

properties (Constant)
    prefix@char scalar = 'P';
end

properties (SetObservable)
    Handheld@logical = true;        % Handhelds have Detector on the top
    SensorType@char = 'msot2';

    ViewMSOTv4Compatible@logical = false;        % be compatible with coordinate definition in ViewMSOT > v4
    
    % geometry
    Curvature@double = 0e-3;               % radius of curvature (in elevation)
    Spacing@double = 0.1e-3;                % Inter-Element Spacing (in m)
    Height@double = 15e-3;                   % Height of Elements

    % casing related properties
    AxialOffset = 0.000;                    % Distance of Membrane from COR
    LightSpotSize = 0.00;                   % Width of LightSpot at Membrane (FWHM, in plane)
    LightSpotHeight = 0.000;                % Height of LightSpot at Membrane (FWHM, out of plane)
    IlluminationCenter = [0 0 0];           % Center of Fiber Bundle
    IlluminationWidth = 0e-2;               % Width of Fiber Bundle
    IlluminationHeight = 0e-3;              % Height of Fiber Bundle
    IlluminationAngle = 0;                  % Inclination Angle of Illumination in degrees
    
    % Simulation Related Properties
    SimulationFoV = [0.06 0.06 0.06]*1e-3;

    % Couplant Properties
    CouplantSOS = 0;                  % Couplant Speed Of Sound
    CouplantCOR = 0;                  % Couplant Center of Rotation
    CouplantCurvature = 0;            % Couplant Radius of Curvature

    % acoustic properties
    CenterFrequency@double = 0e6;     % Center Frequency in Hz (default: 0MHz)
    Bandwidth@double = 0;             % Relative T/R Bandwidth (default: 0.6)
    PhaseShift@double = 0.0;          % Phase Shift (default: 0)
    PhaseSlope@double = 0;            % Phase Slope (default: 0)
    ImpulseResponse@double vector = [];  % Impulse Response (will be computed if not set) 
    
    % Sensitivity maps
    OA_sens_map = {};
    US_sens_map = {};
    
    % New probe file parameters for SW versions >= 4.0 / 2.0
    SOS_Lower_Slider_Limit@int16 = int16(-50);
    SOS_Upper_Slider_Limit@int16 = int16(50);
    Min_Couplant_COR_Offset@double = 0;
    Max_Couplant_COR_Offset@double = 0;
    ADC_Gain@int16 = int16(48);
    ADC_Gain_US@int16 = int16(48);
    Trigger_Delay@int16 = int16(1);
    Trigger_Delay_US@int16 = int16(163);
    OPUS_XDucer_Type@char = 'Imasonic512e270d';
    OPUS_XDucer_App@char = 'General';
    OPUS_Ax_Offset@double = 40;
    OPUS_Master_Gain@int16 = int16(0);
    OPUS_Repeat_Times@int16 = int16(1);
    OPUS_Flip_Image@logical = false;
    OPUS_Acquire@logical = true;
    
end

% Observable. protected Properties - influence the hash code
properties (SetAccess = protected, SetObservable)
    CustomImpulseResponse@logical = false;
end

% Observable Properties - influence the hash code
properties (Abstract = true, SetObservable)
    Sensors@double;           % Cartesian Coordinates of the Sensors (x,y,z)
end

% Don't influence the hash code
properties (Abstract = true)
    Radius@double;            % radius of detector (only Parametric)
    Coverage@double;          % coverage in degree (only Parametric)
    NumDetectors@uint16;      % Number of Elements
    Angles@double;            % Angles (only Parametric)
    StartAngle@double;        % Angle of first Element, read only (only Parametric)
    StepAngle@double;         % Step Angle, read only (only Parametric)
    EndAngle@double;          % End Angle, read only (only Parametric)
    Pitch@double;             % Width of an Element, in Plane (only Parametric)
end

properties (Abstract = true, Constant = true)
    Type@char vector;                   % Type of Probe Definition (Parametric or Coordinates)
    GeometryDefinitionType@char vector; % Type of Probe Definition (Parametric or Coordinates)
    % Two properties coexist to maintain field names of the legacy scheme
end

properties (Constant)
    version@uint8 = uint8(1);
end

properties (Dependent = true)
    N@uint16;                 % Number of Elements
end

methods (Abstract)
    str = writeXML_(obj)
    obj = readXML_(obj,dom)
end

methods
    % *****************************************************************
    %% CONSTRUCTOR
    function obj = msotProbe(ID,DID,varargin)
        % Create MSOT Probe Object
        %
        % msotProbe(ProbeID,DesignID)
        obj@msotSettings(varargin{:});
        obj.ProbeID = ID;
        obj.DesignID = DID;
    end
    
    % *****************************************************************
    %% ACCESSORS
    function N_ = get.N(obj)
        N_ = obj.NumDetectors;
    end
    
    function symm = isSymmetric(obj)
        % returns TRUE if Symmetry of Geometry is within numeric precision (eps)
        S = obj.Sensors;
        symm = ~logical(sum(uint16(diff(abs([S(1:obj.N/2,1) flipud(S(obj.N/2+1:end,1))]),1,2)>eps)) + ...
            sum(uint16(diff(abs([S(1:obj.N/2,2) flipud(S(obj.N/2+1:end,2))]),1,2)>eps)));
    end
    
    function impresp = get.ImpulseResponse(obj)
        if obj.CustomImpulseResponse
            impresp = obj.ImpulseResponse;
        else
            if obj.CenterFrequency == 0 || obj.Bandwidth == 0
                impresp = [];
            else
                IR = msotEIR(obj.CenterFrequency,obj.Bandwidth,obj.PhaseShift);
                impresp = IR.tdom;
            end
        end
    end
    
    function set.ImpulseResponse(obj,I)
        obj.ImpulseResponse = I;
        obj.CustomImpulseResponse = ~isempty(I); 
    end
    
    function set.OA_sens_map(obj, OA_map)
        obj.OA_sens_map = OA_map;
    end
    
    function set.US_sens_map(obj, US_map)
        obj.US_sens_map = US_map;
    end
    
    function writeXML(obj,varargin)
        % write Probe XML File
        %
        % The resulting file contains all information required to
        % characterise the probe and can be read using readXML on an
        % instantiated class, or the static msotProbe.load(FileName) function.
        % 
        % Per Default, an output filename based on the following schema is
        % generated:  <ProbeID>.probe
        %
        % Parameters:
        % * FileName: (optional) Override default FileName
        if isempty(varargin) || isempty(varargin{1})
            filename = [obj.ProbeID '.probe'];
        else
            filename = varargin{1};
        end
            
        
        % open file
        FID = fopen(filename,'w');
        if FID == -1
            error('Cannot write output file');
        end

        % generate header
        fprintf(FID,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
        fprintf(FID,'<DataModelSensorPositions xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n');
        
        % Probe Information
        fprintf(FID,'\t<ProbeID>%s</ProbeID>\n',obj.ProbeID);
        fprintf(FID,'\t<DesignID>%s</DesignID>\n',obj.DesignID);
        hh = 'false';
        if obj.Handheld, hh = 'true'; end
        fprintf(FID,'\t<Handheld>%s</Handheld>\n',hh);
        v4 = 'false';
        if obj.ViewMSOTv4Compatible, v4 = 'true'; end 
        fprintf(FID,'\t<ViewMSOTv4Compatible>%s</ViewMSOTv4Compatible>\n', v4);
        fprintf(FID,'\t<Type>%s</Type>\n',obj.Type);

        % Implementation Class Output;
        fprintf(FID,'%s',writeXML_(obj));
        fprintf(FID,'\t<NumDetectors>%i</NumDetectors>\n',obj.NumDetectors);
        
        % Acoustic Parameters
        fprintf(FID,'\t<CenterFrequency>%.17f</CenterFrequency>\n',obj.CenterFrequency);
        fprintf(FID,'\t<Bandwidth>%.17f</Bandwidth>\n',obj.Bandwidth);
        fprintf(FID,'\t<PhaseShift>%.17f</PhaseShift>\n',obj.PhaseShift);
        

        % Geometric Parameters
        fprintf(FID,'\t<Curvature>%.17f</Curvature>\n',obj.Curvature);
        fprintf(FID,'\t<Spacing>%.17f</Spacing>\n',obj.Spacing);
        fprintf(FID,'\t<Height>%.17f</Height>\n',obj.Height);
        
        % Casing Related Parameters
        fprintf(FID,'\t<AxialOffset>%.17f</AxialOffset>\n',obj.AxialOffset);
        fprintf(FID,'\t<LightSpotSize>%.17f</LightSpotSize>\n',obj.LightSpotSize);
        fprintf(FID,'\t<LightSpotHeight>%.17f</LightSpotHeight>\n',obj.LightSpotHeight);
        fprintf(FID,'\t<IlluminationCenter><x>%.17f</x><y>%.17f</y><z>%.17f</z></IlluminationCenter>\n',obj.IlluminationCenter);
        fprintf(FID,'\t<IlluminationWidth>%.17f</IlluminationWidth>\n',obj.IlluminationWidth);
        fprintf(FID,'\t<IlluminationHeight>%.17f</IlluminationHeight>\n',obj.IlluminationHeight);
        fprintf(FID,'\t<IlluminationAngle>%.17f</IlluminationAngle>\n',obj.IlluminationAngle);

        % Simulation Related Parameters
        fprintf(FID,'\t<SimulationFoV><x>%.17f</x><y>%.17f</y><z>%.17f</z></SimulationFoV>\n',obj.SimulationFoV);
        
        % Couplant Related Parameters
        fprintf(FID,'\t<CouplantSOS>%.17f</CouplantSOS>\n',obj.CouplantSOS);
        fprintf(FID,'\t<CouplantCOR>%.17f</CouplantCOR>\n',obj.CouplantCOR);
        fprintf(FID,'\t<CouplantCurvature>%.17f</CouplantCurvature>\n',obj.CouplantCurvature);

        % Write Detector Positions
        pos_elements = obj.Sensors;
        fprintf(FID,'\t<Sensors>\n');
        for i=1:size(pos_elements,1)
            fprintf(FID,'\t\t<SensorPoint><x>%.17f</x><y>%.17f</y><z>%.17f</z></SensorPoint>\n',pos_elements(i,:));
        end
        fprintf(FID,'\t</Sensors>\n');
        
        % Impulse Response in Base64 encoding
        % to decode: y = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(int8(x)),'double')
        IR = obj.ImpulseResponse;
        if ~isempty(IR)
            IR_base64 = char(org.apache.commons.codec.binary.Base64.encodeBase64(typecast(IR,'uint8')))';
        else
            IR_base64 = '';
        end
        fprintf(FID,'\t<ImpulseResponse>%s</ImpulseResponse>\n',IR_base64);
        cir = 'false';
        if obj.CustomImpulseResponse 
            cir = 'true'; 
        end
        fprintf(FID,'\t<CustomImpulseResponse>%s</CustomImpulseResponse>\n',cir);                

        % write file
        fprintf(FID,'</DataModelSensorPositions>\n');

        % close file
        fclose(FID);

%         % open file
%         FID2 = fopen([filename,'.txt'],'w');
% 
%         % 
%         output_cell2 = {};
%         output_cell2{1} = strcat('Radius: ',num2str(Radius,'%.17f'));
%         output_cell2{2} = strcat('Step Angle: ',num2str(StepAngle,'%.17f'));
%         output_cell2{3} = strcat('Start Angle: ',num2str(StartAngle,'%.17f'));
% 
%         % write file
%         fprintf(FID2, '%s\n', output_cell2{:});
% 
%         % close file
%         fclose(FID2);
    end

    function writeXML_v2(obj,varargin)
        % write Probe XML File for vMc v2.0
        %
        % The resulting file contains all information required to
        % characterise the probe and can be read using readXML on an
        % instantiated class, or the static msotProbe.load(FileName) function.
        % 
        % Per Default, an output filename based on the following schema is
        % generated:  <ProbeID>.probe
        %
        % Parameters:
        % * FileName: (optional) Override default FileName
        if isempty(varargin) || isempty(varargin{1})
            filename = [obj.ProbeID '.probe'];
        else
            filename = varargin{1};
        end
            
        
        % open file
        FID = fopen(filename,'w');
        if FID == -1
            error('Cannot write output file');
        end

        % generate header
        fprintf(FID,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
        fprintf(FID,'<DataModelProbeRedistributable xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n');
        
        % Probe Information
        fprintf(FID,'\t<PType>%s</PType>\n',obj.ProbeID);
                
        hh = 'false';
        if obj.Handheld 
            hh = 'true'; 
        end
        fprintf(FID,'\t<Handheld>%s</Handheld>\n',hh);
        
        fprintf(FID,'\t<SensorType>%s</SensorType>\n',obj.SensorType);
        fprintf(FID,'\t<GeometryDefinitionType>%s</GeometryDefinitionType>\n',obj.GeometryDefinitionType);

        % Implementation Class Output;
        fprintf(FID,'%s',writeXML_(obj)); % Writing Radius and Coverage!! --> Abstract method implemented in msotProbeCoordinates/Parametric
        fprintf(FID,'\t<NumDetectors>%i</NumDetectors>\n',obj.NumDetectors);
        
        % Acoustic Parameters
        fprintf(FID,'\t<CenterFrequency>%.17f</CenterFrequency>\n',obj.CenterFrequency);
        fprintf(FID,'\t<Bandwidth>%.17f</Bandwidth>\n',obj.Bandwidth);
        fprintf(FID,'\t<PhaseShift>%.17f</PhaseShift>\n',obj.PhaseShift);
        
        % Geometric Parameters
        fprintf(FID,'\t<Curvature>%.17f</Curvature>\n',obj.Curvature);
        fprintf(FID,'\t<Spacing>%.17f</Spacing>\n',obj.Spacing);
        fprintf(FID,'\t<Height>%.17f</Height>\n',obj.Height);
        
        % Casing Related Parameters
        fprintf(FID,'\t<AxialOffset>%.17f</AxialOffset>\n',obj.AxialOffset);
        fprintf(FID,'\t<LightSpotSize>%.17f</LightSpotSize>\n',obj.LightSpotSize);
        fprintf(FID,'\t<LightSpotHeight>%.17f</LightSpotHeight>\n',obj.LightSpotHeight);
        fprintf(FID,'\t<IlluminationCenter><x>%.17f</x><y>%.17f</y><z>%.17f</z></IlluminationCenter>\n',obj.IlluminationCenter);
        fprintf(FID,'\t<IlluminationWidth>%.17f</IlluminationWidth>\n',obj.IlluminationWidth);
        fprintf(FID,'\t<IlluminationHeight>%.17f</IlluminationHeight>\n',obj.IlluminationHeight);
        fprintf(FID,'\t<IlluminationAngle>%.17f</IlluminationAngle>\n',obj.IlluminationAngle);

        % Simulation Related Parameters
        fprintf(FID,'\t<SimulationFoV><x>%.17f</x><y>%.17f</y><z>%.17f</z></SimulationFoV>\n',obj.SimulationFoV);
        
        % Couplant Related Parameters
        fprintf(FID,'\t<CouplantSOS>%.17f</CouplantSOS>\n',obj.CouplantSOS);
        fprintf(FID,'\t<CouplantCOR>%.17f</CouplantCOR>\n',obj.CouplantCOR);
        fprintf(FID,'\t<CouplantCurvature>%.17f</CouplantCurvature>\n',obj.CouplantCurvature);

        % Write Detector Positions
        pos_elements = obj.Sensors;
        fprintf(FID,'\t<Sensors>\n');
        for i=1:size(pos_elements,1)
            fprintf(FID,'\t\t<SensorPoint><x>%.17f</x><y>%.17f</y><z>%.17f</z></SensorPoint>\n',pos_elements(i,:));
        end
        fprintf(FID,'\t</Sensors>\n');
        
        % Impulse Response in Base64 encoding
        % to decode: y = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(int8(x)),'double')
        IR = obj.ImpulseResponse;
        if ~isempty(IR)
            IR_base64 = char(org.apache.commons.codec.binary.Base64.encodeBase64(typecast(IR,'uint8')))';
        else
            IR_base64 = '';
        end
        fprintf(FID,'\t<ImpulseResponse>%s</ImpulseResponse>\n',IR_base64);
        
        cir = 'false';
        if obj.CustomImpulseResponse 
            cir = 'true'; 
        end
        fprintf(FID,'\t<CustomImpulseResponse>%s</CustomImpulseResponse>\n',cir);                

        
        
        % OA sensitivity map in Base64 encoding
        % to decode: y = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(int8(x)),'double')
        OA = obj.OA_sens_map;
        if ~isempty(OA)
            % First step before encoding is to create a vector of uint8
            % representation (stored as uint8 itself!) of the mixed type content of the sens file.
            % The different entries are gathered in the same vector but the
            % length of their byte representation is preserved
            
            % To pre-allocate the vector of uint8, compute the number of
            % bytes used for representation of the content of the .sens
            % file
            byte_size = 0;
            for i = 1:length(OA)
                entry = OA{i};
                data = whos('entry');
                byte_size = byte_size + data.bytes;
                clear entry data
            end
            vect_pre_coding = uint8(zeros(1, byte_size));
            
            bytes_converted = 0;
            for i = 1:length(OA)
                uint8_representation = typecast(OA{i}, 'uint8');
                vect_pre_coding(1+bytes_converted : bytes_converted+length(uint8_representation)) = uint8_representation;
                bytes_converted = bytes_converted + length(uint8_representation);
            end
            
            % Encoding can now take place
            OA_base64 = char(org.apache.commons.codec.binary.Base64.encodeBase64(vect_pre_coding))';
        else
            OA_base64 = '';
        end
        fprintf(FID,'\t<OptoacousticSensitivityMap>%s</OptoacousticSensitivityMap>\n',OA_base64);
        
        
        
        % US sensitivity map in Base64 encoding
        % to decode: y = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(int8(x)),'double')
        US = obj.US_sens_map;
        if ~isempty(US)
            % First step before encoding is to create a vector of uint8
            % representation (stored as uint8 itself!) of the mixed type content of the sens file.
            % The different entries are gathered in the same vector but the
            % length of their byte representation is preserved
            
            % To pre-allocate the vector of uint8, compute the number of
            % bytes used for representation of the content of the .sens
            % file
            byte_size = 0;
            for i = 1:length(US)
                entry = US{i};
                data = whos('entry');
                byte_size = byte_size + data.bytes;
                clear entry data
            end
            vect_pre_coding = uint8(zeros(1, byte_size));
            
            bytes_converted = 0;
            for i = 1:length(US)
                uint8_representation = typecast(US{i}, 'uint8');
                vect_pre_coding(1+bytes_converted : bytes_converted+length(uint8_representation)) = uint8_representation;
                bytes_converted = bytes_converted + length(uint8_representation);
            end
            
            % Encoding can now take place
            US_base64 = char(org.apache.commons.codec.binary.Base64.encodeBase64(vect_pre_coding))';
        else
            US_base64 = '';
        end
        fprintf(FID,'\t<UltrasoundSensitivityMap>%s</UltrasoundSensitivityMap>\n',US_base64);
        
        % New probe file parameters for SW versions >= 4.0 / 2.0
        fprintf(FID,'\t<SpeedOfSoundLowerSliderLimit>%i</SpeedOfSoundLowerSliderLimit>\n',obj.SOS_Lower_Slider_Limit);
        fprintf(FID,'\t<SpeedOfSoundUpperSliderLimit>%i</SpeedOfSoundUpperSliderLimit>\n',obj.SOS_Upper_Slider_Limit);
        fprintf(FID,'\t<MinCouplantCenterOfRotationOffset>%.17f</MinCouplantCenterOfRotationOffset>\n',obj.Min_Couplant_COR_Offset);
        fprintf(FID,'\t<MaxCouplantCenterOfRotationOffset>%.17f</MaxCouplantCenterOfRotationOffset>\n',obj.Max_Couplant_COR_Offset);
        fprintf(FID,'\t<ADCGain>%i</ADCGain>\n',obj.ADC_Gain);
        fprintf(FID,'\t<ADCGainUS>%i</ADCGainUS>\n',obj.ADC_Gain_US);
        fprintf(FID,'\t<TriggerDelay>%i</TriggerDelay>\n',obj.Trigger_Delay);
        fprintf(FID,'\t<TriggerDelayUS>%i</TriggerDelayUS>\n',obj.Trigger_Delay_US);
        fprintf(FID,'\t<OpusXDucerType>%s</OpusXDucerType>\n',obj.OPUS_XDucer_Type);
        fprintf(FID,'\t<OpusXDucerApp>%s</OpusXDucerApp>\n',obj.OPUS_XDucer_App);
        fprintf(FID,'\t<OpusAxOffset>%.17f</OpusAxOffset>\n',obj.OPUS_Ax_Offset);
        fprintf(FID,'\t<OpusMasterGain>%i</OpusMasterGain>\n',obj.OPUS_Master_Gain);
        fprintf(FID,'\t<OpusRepeatTimes>%i</OpusRepeatTimes>\n',obj.OPUS_Repeat_Times);
        
        OPUS_FI = 'false';
        if obj.OPUS_Flip_Image
            OPUS_FI = 'true'; 
        end
        fprintf(FID,'\t<OpusFlipImage>%s</OpusFlipImage>\n', OPUS_FI);
        
        OPUS_A = 'false';
        if obj.OPUS_Acquire
            OPUS_A = 'true'; 
        end
        fprintf(FID,'\t<OpusAcquire>%s</OpusAcquire>\n', OPUS_A);
        
        % write file
        fprintf(FID,'</DataModelProbeRedistributable>\n');

        % close file
        fclose(FID);
    end

    function writeIRFFile(obj,varargin)
        % write file
        if isempty(varargin) || isempty(varargin{1})
            filename = [obj.ProbeID '.irf'];
        else
            filename = varargin{1};
        end
        
        FID = fopen(filename,'w');
        fwrite(FID,obj.ImpulseResponse,'double');
        fclose(FID);  
    end
    
    function [schx, schz] = ProdigyINI(obj,varargin)
        % plot the Detector Arrangement for the Prodigy Configuration
        %
        % copy the function output to the [RUCT] secrion of the Prodigy INI
        % File.
        % NOTE: This is also contained in the Probe XML File written by
        % writeXML
        if isempty(varargin) || isempty(varargin{1})
            filename = [obj.ProbeID '_Prodigy.txt'];
        else
            filename = varargin{1};
        end
        
        % open file
        FID = fopen(filename,'w');
        if FID == -1
            error('Cannot write output file');
        end
        
        if  ~obj.ViewMSOTv4Compatible && ~obj.Handheld % pre-clinical 1st ch +X +Z
            pos_elements_us = obj.Sensors; 
            pos_elements_us(:,2) = (-1)*pos_elements_us(:,2);
        elseif obj.ViewMSOTv4Compatible && ~obj.Handheld % pre-clinical 1st ch +X +Z
            pos_elements_us = [(-1)*obj.Sensors(:,1) obj.Sensors(:,3) obj.Sensors(:,2) ];
        elseif ~obj.ViewMSOTv4Compatible && obj.Handheld % clinical 1st ch
            pos_elements_us = [obj.Sensors(:,1) (-1)*obj.Sensors(:,2) obj.Sensors(:,3) ];
        elseif obj.ViewMSOTv4Compatible && obj.Handheld % clinical 1st ch
            pos_elements_us = [(-1)*obj.Sensors(:,1) (-1)*obj.Sensors(:,3) obj.Sensors(:,2) ];    
        end
        digits = 17;
        schx='[';
        schz='[';
        for ii=1:obj.NumDetectors-1
            schx = strcat(schx,num2str(pos_elements_us(ii,1),digits),',');
            schz = strcat(schz,num2str(pos_elements_us(ii,2),digits),',');  
        end
            schx = strcat('Chx=0,1,',num2str(obj.NumDetectors),',1:',schx,num2str(pos_elements_us(obj.NumDetectors,1),digits),']');
            schz = strcat('Chz=0,1,',num2str(obj.NumDetectors),',1:',schz,num2str(pos_elements_us(obj.NumDetectors,2),digits),']');             

        % copy to [RUCT] field of ini file
        disp(schx);
        disp(schz);
        
        % write file
        fprintf(FID,'++ Copy the below to [RUCT] section of ini file ++\r\n');
        fprintf(FID,'%s\r\n',schx);
        fprintf(FID,'%s\r\n',schz);

        % close file
        fclose(FID);

    end
    
    function [fig, ax1] = plot(obj,varargin)
        % plot the Detector Arrangement
        %
        % Optional Parameters:
        % - Axis Handles
        
        % Sanity checks
        if isempty(obj.Sensors), error('No Sensor Defined'); end
        if std(obj.Sensors(:,3)) > eps
            error('Cannot plot three-dimensional Detector (yet)');
        end
        
        % Constants (maybe make configurable at some point?
        par.TransducerColor = [     0    0.4470    0.7410];
        par.MembraneColor = [0.4940    0.1840    0.5560];
        par.ScaleColor = [0.5 0.5 0.5];
        par.LightSpotColor = [0.8500    0.1250    0.0980];
        par.PitchColor = [    0.4660    0.6740    0.1880];
        par.SimulationColor = [0.5 0.5 0.5];
        
        % Figure Creation
        if nargin >= 2 && isvalid(varargin{1})
            ax1 = varargin{1};
            axes(ax1);      
            fig = gcf;      % just get figure handle as well
            if nargin >= 3
                ax2 = varargin{2};
            else
                ax2 = [];
            end
        else
            fig = figure('Position',[100 100 1200 800]);
            ax1 = axes('Units','pixels','Position',[40 40 700 700]);
            if obj.Height > 0
                ax2 = axes('Units','pixels','Position',[780 40 round(700/3) 700]);
            else
                ax2 = [];
            end
        end
        
        % Center of Rotation
        axes(ax1);
        plot(0,0,'x','Color','k','MarkerSize',6); hold on;
        text(0,-2,'(0,0)','Color','k','FontSize',8,'VerticalAlign','top','HorizontalAlign','center');
        
        % Sensors (in mm!)
        Sx = obj.Sensors(:,1)*1e3;
        Sy = -obj.Sensors(:,2)*1e3;
        ang_ = atan2(Sy,Sx);
        rad_ = median(sqrt(Sx.^2 + Sy.^2));

        % Axis Dimension
        xl = [-rad_*1.5 rad_*1.5];
        xlim(xl);
        yl = [-rad_*1.5 rad_*1.5];
        ylim(yl);
        set(ax1,'XTick',[]);
        set(ax1,'YTick',[]);
        
        % Probe ID
%         text(xl(2)*0.9,yl(1)*0.9,{sprintf('\\hfill\\textbf{ProbeID: %s}',obj.ProbeID),sprintf('DesignID: %s',obj.DesignID)},'HorizontalAlign','right','VerticalAlign','bottom','FontSize',14,'Interpreter','latex');
        text(xl(2)*0.9,yl(1)*0.9,{sprintf('\\bfProbeID: %s',obj.ProbeID),sprintf('\\rmDesignID: %s',obj.DesignID)},'HorizontalAlign','right','VerticalAlign','bottom','FontSize',14);
        
        % membrane (TODO: Use Couplant Geometry)
        if obj.AxialOffset > 0
            line(xlim(),[1 1]*obj.AxialOffset*1e3,'Color',par.MembraneColor,'LineStyle','--');
            text(min(Sx),obj.AxialOffset*1e3,'Membrane','Color',par.MembraneColor,'HorizontalAlign','right','VerticalAlign','bottom');
            line([1 1]*min(Sx),[0 obj.AxialOffset]*1e3,'Color',par.MembraneColor,'LineStyle','-');
            line([0.95 1.05]*1*min(Sx),[0 0],'Color',par.MembraneColor,'LineStyle','-');
            text(1*min(Sx),obj.AxialOffset*0.5e3,sprintf('%.1fmm  ',obj.AxialOffset*1e3),'Color',par.MembraneColor,'LineStyle','-','VerticalAlign','middle','HorizontalAlign','right');

            if obj.LightSpotSize > 0
                line([-0.5 0.5]*obj.LightSpotSize*1e3,[1 1]*obj.AxialOffset*1e3,'Color',par.LightSpotColor,'LineWidth',3);
                text(0,obj.AxialOffset*1e3,sprintf('FWHM = %1.fmm',obj.LightSpotSize*1e3),'Color',par.LightSpotColor,'VerticalAlign','top','HorizontalAlign','center');
            end
        end
        
        % fiber bundle
        if sum(obj.IlluminationCenter) > 0 && obj.IlluminationWidth > 0
           line(obj.IlluminationCenter(1)*1e3+[-0.5 0.5]*obj.IlluminationWidth*1e3,[1 1]*obj.IlluminationCenter(2)*1e3,'LineWidth',5,'Color',par.LightSpotColor); 
           text(0,obj.IlluminationCenter(2)*1e3+1,sprintf('%1.fmm',obj.IlluminationWidth*1e3),'Color',par.LightSpotColor,'VerticalAlign','bottom','HorizontalAlign','center');
        end

        % plot Detector
        Sx_ = flip(rad_.*1.25.*cos(ang_));
        Sy_ = flip(rad_.*1.25.*sin(ang_));
        plot([Sx;Sx_;Sx(1)],[Sy;Sy_;Sy(1)],'Color',par.TransducerColor);hold on;
        plot(Sx,Sy,'.','Color',par.TransducerColor,'MarkerSize',5);
        plot([Sx(1) 0 Sx(end)],[Sy(1) 0 Sy(end)],'LineStyle',':','Color',par.TransducerColor);
        text(Sx(end)/2,Sy(end)/2,sprintf('r = %.1fmm',obj.Radius*1e3),'Color',par.TransducerColor,'VerticalAlign','middle','HorizontalAlign','center','Rotation',atan2(Sy(end),Sx(end))/pi*180,'BackgroundColor','white');
        text(0,(-1 + 2*obj.Handheld)*rad_/3,{sprintf('\\bf%.1f°',obj.Coverage),sprintf('\\rm%i Elements',obj.N)},'HorizontalAlign','center','VerticalAlign','middle','Color',par.TransducerColor,'FontSize',12);
%         text((rad_*0.8)*cos(median(ang_)),(rad_*0.8)*sin(median(ang_)),sprintf('h = %.1fmm\nr_c = %.1fmm',obj.Height*1e3,obj.Curvature*1e3),...
%             'HorizontalAlign','center',...
%             'VerticalAlign','middle','Color',par.TransducerColor);
%         text(Sx(1)*0.9,Sy(1),sprintf('Pitch = %.3fmm',obj.Pitch*1e3),'Color',par.TransducerColor,'VerticalAlign','middle')
        text(Sx(1)+1,Sy(1)-1,'1','Color',par.TransducerColor,'VerticalAlign','top','HorizontalAlign','left','BackgroundColor','white','Margin',1);
        text(Sx(end)-1,Sy(end)-1,num2str(obj.N),'Color',par.TransducerColor,'VerticalAlign','top','HorizontalAlign','right','BackgroundColor','white','Margin',1);

        % Pitch
        [~,pel] = min(abs(Sx));
        plot(Sx(pel:pel+1),Sy(pel:pel+1),'Color',par.PitchColor,'LineWidth',2);
        text(Sx(pel),Sy(pel)*0.95,sprintf('Pitch = %.3fmm',obj.Pitch*1e3),'Color',par.PitchColor,'VerticalAlign','middle','Rotation',double(~obj.Handheld)*180-90);
        
        % Transducer scale
        line([min(Sx) max(Sx)],rad_*1.4*[1 1],'Color',par.ScaleColor);
        line(min(Sx)*[1 1],rad_*1.4*[0.98 1.02],'Color',par.ScaleColor);
        line(max(Sx)*[1 1],rad_*1.4*[0.98 1.02],'Color',par.ScaleColor);
        text(0,rad_*1.4,sprintf('%.1fmm',diff([min(Sx),max(Sx)])),'Color',par.ScaleColor,'FontSize',10,'VerticalAlign','middle','HorizontalAlign','center','BackgroundColor','white');
        
        % Simulation FoV
        if sum(obj.SimulationFoV) > 0
            plot([-1 -1 1 1 -1]*obj.SimulationFoV(1)*0.5e3,[1 -1 -1 1 1]*obj.SimulationFoV(2)*0.5e3,'Color',par.SimulationColor,'LineStyle','--');
            text(-obj.SimulationFoV(1)*0.5e3,-obj.SimulationFoV(1)*0.5e3,sprintf('Simulation Field of View (%ix%ix%imm)',round(obj.SimulationFoV*1e3)),'Color',par.SimulationColor,'VerticalAlign','top');
        end
        
        % coordinate system
        line(xl(1)*0.9+[0 10],yl(1)*0.9+[0 0],'Color','k'); 
        line(xl(1)*0.9+[0 0],yl(1)*0.9+[0 10],'Color','k'); 
        text(0.9*xl(1),0.7*yl(1),' y','FontSize',8,'HorizontalAlign','left','VerticalAlign','top');
        text(xl(1)*0.9+10,0.9*yl(1),'x','FontSize',8,'HorizontalAlign','right','VerticalAlign','top');
        text(xl(1)*0.9+5,0.9*yl(1),'1cm','FontSize',8,'HorizontalAlign','center','VerticalAlign','bottom');
        
        % second axis
        if ~isempty(ax2)
            axes(ax2);
            plot(0,0,'x','Color','k','MarkerSize',6); hold on;
            text(0,0,'(0,0)','Color','k','FontSize',8,'VerticalAlign','top');
            
            % Axis Dimension
            xlim(xl/3);
            yl = [-rad_*1.5 rad_*1.5];
            ylim(yl);
            set(ax2,'XTick',[]);
            set(ax2,'YTick',[]);
            
            rad_c = obj.Curvature*1e3;
            rdiff = rad_ - rad_c;
            ang_c = asin(obj.Height/2/obj.Curvature);
            ang_c = linspace(-ang_c,ang_c,100);
            Sy_c = rdiff+rad_c*cos(ang_c);
            if ~obj.Handheld, Sy_c = -Sy_c; end
            Sz_c = rad_c*sin(ang_c);
            plot(Sz_c,Sy_c,'.','Color',par.TransducerColor,'MarkerSize',5);
            plot([Sz_c(1) 0 Sz_c(end)],[Sy_c(1) rdiff Sy_c(end)],...
                'LineStyle',':','Color',par.TransducerColor);
            text(Sz_c(end)/2,Sy_c(end)/2,sprintf('r_c = %.1fmm',obj.Curvature*1e3),'Color',par.TransducerColor,'VerticalAlign','middle','HorizontalAlign','center','Rotation',90-atan2(Sz_c(end),Sy_c(end))/pi*180,'BackgroundColor','white');
            
            % membrane (TODO: Use Couplant Geometry)
            if obj.AxialOffset > 0
                line(ylim(),[1 1]*obj.AxialOffset*1e3,'Color',par.MembraneColor,'LineStyle','--');
                text(max(Sz_c)+10,obj.AxialOffset*1e3,' Membrane','Color',par.MembraneColor,'HorizontalAlign','right','VerticalAlign','bottom');
    %             line([1 1]*min(Sx),[0 obj.AxialOffset]*1e3,'Color',par.MembraneColor,'LineStyle','-');
    %             line([0.95 1.05]*1*min(Sx),[0 0],'Color',par.MembraneColor,'LineStyle','-');
    %             text(1*min(Sx),obj.AxialOffset*0.5e3,sprintf('%.1fmm  ',obj.AxialOffset*1e3),'Color',par.MembraneColor,'LineStyle','-','VerticalAlign','middle','HorizontalAlign','right');
                if obj.LightSpotHeight > 0
                    line([-0.5 0.5]*obj.LightSpotHeight*1e3,[1 1]*obj.AxialOffset*1e3,'Color',par.LightSpotColor,'LineWidth',3);
                    text(0,obj.AxialOffset*1e3,sprintf('FWHM = %1.fmm',obj.LightSpotHeight*1e3),'Color',par.LightSpotColor,'VerticalAlign','top','HorizontalAlign','center');
                end            
            end
            
            % fiber bundle
            if sum(obj.IlluminationCenter) > 0 && obj.IlluminationWidth > 0
               ang_f = -obj.IlluminationAngle/180*pi;
               Sz_f = obj.IlluminationCenter(3)*1e3+obj.IlluminationHeight*1e3*0.5*[-1,1]*cos(ang_f);
               Sy_f = obj.IlluminationCenter(2)*1e3+obj.IlluminationHeight*1e3*0.5*[-1,1]*sin(ang_f);
               line(Sz_f,Sy_f,'LineWidth',5,'Color',par.LightSpotColor); 
               text(obj.IlluminationCenter(3)*1e3,obj.IlluminationCenter(2)*1e3-1,sprintf('%1.fmm',obj.IlluminationHeight*1e3),'Color',par.LightSpotColor,'VerticalAlign','top','HorizontalAlign','center','Rotation',-obj.IlluminationAngle);
            end

            % Transducer scale
            line([min(Sz_c) max(Sz_c)],rdiff+rad_c*1.4*[1 1],'Color',par.ScaleColor);
            line(min(Sz_c)*[1 1],rad_*1.4*[0.98 1.02],'Color',par.ScaleColor);
            line(max(Sz_c)*[1 1],rad_*1.4*[0.98 1.02],'Color',par.ScaleColor);
            text(0,rad_*1.4,sprintf('%.1fmm',diff([min(Sz_c),max(Sz_c)])),'Color',par.ScaleColor,'FontSize',10,'VerticalAlign','middle','HorizontalAlign','center','BackgroundColor','white');%,'Rotation',90);

            % Simulation FoV
            if sum(obj.SimulationFoV) > 0
                plot([-1 -1 1 1 -1]*obj.SimulationFoV(3)*0.5e3,[1 -1 -1 1 1]*obj.SimulationFoV(2)*0.5e3,'Color',par.SimulationColor,'LineStyle','--');
%                 text(-obj.SimulationFoV(3)*0.5e3,-obj.SimulationFoV(1)*0.5e3,sprintf('Simulation Field of View (%ix%ix%imm)',round(obj.SimulationFoV*1e3)),'Color',par.SimulationColor,'VerticalAlign','top');
            end

            % coordinate system
            line(xl(1)/3+5+[0 10],yl(1)*0.9+[0 0],'Color','k'); 
            line(xl(1)/3+5+[0 0],yl(1)*0.9+[0 10],'Color','k'); 
            text(0.9*xl(1)/3,0.7*yl(1),' y','FontSize',8,'HorizontalAlign','left','VerticalAlign','top');
            text(xl(1)/3+15,0.9*yl(1),'z','FontSize',8,'HorizontalAlign','right','VerticalAlign','top');
            text(xl(1)/3+10,0.9*yl(1),'1cm','FontSize',8,'HorizontalAlign','center','VerticalAlign','bottom');
        end
    end
end


%% ************************************************************************
% Listener and Property Modification

methods (Access = protected)
    function createListeners(obj)
        % Associate all events with PostSet listeners to re-calculate
        % hash and notify potential Listeners
        mc = metaclass(obj);
        proplist = mc.PropertyList;
        ind = find(~ismember({proplist.Name},{'version','hashCode','ID','sID','logger','creation','undoProp','undoValue'}));
        for ix = ind
            % Ignore constants and not observable properties (will fail...)
            if proplist(ix).Constant || ~proplist(ix).SetObservable
                continue 
            end 
            addlistener(obj, proplist(ix).Name, 'PreSet', @obj.propertyNotYetModified);
            addlistener(obj, proplist(ix).Name, 'PostSet', @obj.propertyModified);
        end
    end
end



methods (Static = true)
    function obj = readXML(varargin)
        % Read any type of Probe XML file
        % Will return the correct msotProbe object, based on the type
        % specified in the XML "Type" or "GeometryDefinitionType" Node. 
        %
        % obj = msotProbe.readXML(filename[,obj])
        %   read XML File given in filename, and instantiate Probe object
        %   or populate Probe Object in obj
        %
        % obj = msotProbe.readXML(dom)
        %   assign properties from DOM object dom to a new msotProbe object
        %   of correct type or populate msotProbe object in obj
        
        if numel(varargin) == 0
            error('Please specify at least a filename or DOM object');
        end
        
        % Load or use DOM object
        if isa(varargin{1},'char')
            dom = msotProbe.loadXML(varargin{1}); %dom = patdoc.getDataModelSensorPositions(); OR patdoc.getDataModelProbeRedistributable();
            % Where patdoc is the javaMethod from msotbeans DataModelSensorPositions OR DataModelProbeRedistributable
        elseif isa(varargin{1},'DOM')
            dom = varargin{1};
        end
        
        % Determine the type of .probe file being loaded (either _v2 or old type)
        object_class = class(dom);
        if strcmp(object_class(34:57), 'DataModelSensorPositions')
            file_type = 'old';
        else
            file_type = 'new';
        end
        
        % Read Minimum Properties from DOM (possibly required to instantiate a fresh probe object)
        if strcmp(file_type, 'old')
            prid = char(dom.getProbeID());
            did = char(dom.getDesignID());
            type = char(dom.getType());
        else
            prid = char(dom.getPType());
            did = 'unkownDesign';
            type = char(dom.getGeometryDefinitionType());
        end
                
        s = dom.getSensors.getSensorPointArray;
        S = zeros(numel(s),3);
        for jj = 1:numel(s)
            S(jj,:) = double([s(jj).getX s(jj).getY s(jj).getZ]);
        end
        
        HH = dom.getHandheld;
        v4 = logical(0);
        if strcmp(file_type, 'old')
            if dom.isSetViewMSOTv4Compatible 
                v4 = dom.getViewMSOTv4Compatible; 
            end
        else
            v4 = true;
        end
        
        % Update essential class arguments
        if nargin == 2 && isa(varargin{2}, 'msotProbe')
            obj = varargin{2};
            if strcmp(file_type, 'old') && ~strcmp(obj.Type, type)
                error('Detector Type Mismatch');
            elseif strcmp(file_type, 'new') && ~strcmp(obj.GeometryDefinitionType, type)
                error('Detector GeometryDefinitionType Mismatch');
            end
            
            obj.ProbeID = prid;
            obj.DesignID = did;
            % viewMSOTv4Compatible and Handheld fields are added later
            
            if isa(obj,'msotProbeParametric')
                obj.Radius = double(dom.getRadius);
                obj.Coverage = double(dom.getCoverage);
                obj.NumDetectors = uint16(dom.getNumDetectors);
                % Useless to write obj.Sensors = S since the positions are calculated from R, C and N
                
            elseif isa(obj,'msotProbeCoordinates')
                obj.Sensors = S;
                % Useless to specify R, C and N since those are calculated from the sensor positions
            end
            
        else
        % Or instantiate new object    
            switch type
                case 'Parametric'
                    R = double(dom.getRadius);
                    C = double(dom.getCoverage);
                    N = uint16(dom.getNumDetectors);
                    obj = msotProbeParametric(prid, did, R, C, N, S, HH, v4);
                
                % In all other cases use Coordinates - doesn't need anything else in terms of attributes
                otherwise
                    obj = msotProbeCoordinates(prid, did, S, v4);
            end
        end
        
        % --> Write in obj all the parameters read in the probe file
        obj.Handheld = HH;
        obj.ViewMSOTv4Compatible = v4;
        
        % Impulse Reponse
        obj.CustomImpulseResponse = dom.getCustomImpulseResponse; % dom.getField returns a logical 0 or 1 for fields of type boolean
        if obj.CustomImpulseResponse
            IR = char(dom.xgetImpulseResponse.getStringValue);
            % Convert IR back from base64 encoding to 1D signal
            obj.ImpulseResponse = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(uint8(IR)), 'double');
        else
            obj.ImpulseResponse = [];
        end
        
        % Acoustic
        if dom.isSetCenterFrequency, obj.CenterFrequency = double(dom.getCenterFrequency); end
        if dom.isSetBandwidth, obj.Bandwidth = double(dom.getBandwidth); end
        if dom.isSetPhaseShift, obj.PhaseShift = double(dom.getPhaseShift); end

        % Geometric
        if dom.isSetCurvature, obj.Curvature = double(dom.getCurvature); end
        if dom.isSetSpacing, obj.Spacing = double(dom.getSpacing); end
        if dom.isSetHeight, obj.Height = double(dom.getHeight); end
        
        % Casing Related
        if dom.isSetAxialOffset, obj.AxialOffset = double(dom.getAxialOffset); end
        if dom.isSetLightSpotSize, obj.LightSpotSize = double(dom.getLightSpotSize); end
        if dom.isSetLightSpotHeight, obj.LightSpotHeight = double(dom.getLightSpotHeight); end
        if dom.isSetIlluminationCenter, obj.IlluminationCenter = double([dom.getIlluminationCenter.getX dom.getIlluminationCenter.getY dom.getIlluminationCenter.getZ])'; end
        if dom.isSetIlluminationWidth, obj.IlluminationWidth = double(dom.getIlluminationWidth); end
        if dom.isSetIlluminationHeight, obj.IlluminationHeight = double(dom.getIlluminationHeight); end
        if dom.isSetIlluminationAngle, obj.IlluminationAngle = double(dom.getIlluminationAngle); end
        
        % Simulation
        if dom.isSetSimulationFoV, obj.SimulationFoV = double([dom.getSimulationFoV.getX dom.getSimulationFoV.getY dom.getSimulationFoV.getZ ])'; end

        % Couplant Related
        if dom.isSetCouplantSOS, obj.CouplantSOS = double(dom.getCouplantSOS); end
        if dom.isSetCouplantCOR, obj.CouplantCOR = double(dom.getCouplantCOR); end
        if dom.isSetCouplantCurvature, obj.CouplantCurvature = double(dom.getCouplantCurvature); end       
        
        % --> Complete reading of the probe file with fields specific to the new format
        if strcmp(file_type, 'new')
            if dom.isSetSensorType, obj.SensorType = char(dom.getSensorType); end
            if dom.isSetSpeedOfSoundLowerSliderLimit, obj.SOS_Lower_Slider_Limit = int16(dom.getSpeedOfSoundLowerSliderLimit); end
            if dom.isSetSpeedOfSoundUpperSliderLimit, obj.SOS_Upper_Slider_Limit = int16(dom.getSpeedOfSoundUpperSliderLimit); end
            if dom.isSetMinCouplantCenterOfRotationOffset, obj.Min_Couplant_COR_Offset = double(dom.getMinCouplantCenterOfRotationOffset); end
            if dom.isSetMaxCouplantCenterOfRotationOffset, obj.Max_Couplant_COR_Offset = double(dom.getMaxCouplantCenterOfRotationOffset); end
            if dom.isSetADCGain, obj.ADC_Gain = int16(dom.getADCGain); end
            if dom.isSetADCGainUS, obj.ADC_Gain_US = int16(dom.getADCGainUS); end
            if dom.isSetTriggerDelay, obj.Trigger_Delay = int16(dom.getTriggerDelay); end
            if dom.isSetTriggerDelayUS, obj.Trigger_Delay_US = int16(dom.getTriggerDelayUS); end
            if dom.isSetOpusXDucerType, obj.OPUS_XDucer_Type = char(dom.getOpusXDucerType); end
            if dom.isSetOpusXDucerApp, obj.OPUS_XDucer_App = char(dom.getOpusXDucerApp); end
            if dom.isSetOpusAxOffset, obj.OPUS_Ax_Offset = double(dom.getOpusAxOffset); end
            if dom.isSetOpusMasterGain, obj.OPUS_Master_Gain = int16(dom.getOpusMasterGain); end
            if dom.isSetOpusRepeatTimes, obj.OPUS_Repeat_Times = int16(dom.getOpusRepeatTimes); end
            if dom.isSetOpusFlipImage, obj.OPUS_Flip_Image = dom.getOpusFlipImage; end
            if dom.isSetOpusAcquire, obj.OPUS_Acquire = dom.getOpusAcquire; end
            
            % Decoding the OA sens map
            if dom.isSetOptoacousticSensitivityMap
                % Because the sensitivity map contains two different data types (both
                % uint16 and single), the base64 decoding has to be performed with care. 
                % Indeed each entry of the file should be decoded according to its own type.
                % The content of the .sens file is stored in a cell array,
                % allowing for mixed data-type (unlike vectors or arrays) with one element per cell.

                code = char(dom.xgetOptoacousticSensitivityMap.getStringValue);
                
                if isempty(code)
                    obj.OA_sens_map = {};
                    % Covers the case when the sens map field exists but is empty
                else
                    decode = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(uint8(code)), 'uint8');
                    sens_map = {};

                    % Decode header
                    sens_map{1} = typecast(decode(1:2), 'uint16')';
                    sens_map{2} = typecast(decode(3:4), 'uint16')';
                    sens_map{3} = typecast(decode(5:6), 'uint16')';
                    sens_map{4} = typecast(decode(7:8), 'uint16')';

                    sens_map{5} = typecast(decode(9:12), 'single')';
                    sens_map{6} = typecast(decode(13:16), 'single')';
                    sens_map{7} = typecast(decode(17:20), 'single')';

                    sens_map{8} = typecast(decode(21:22), 'uint16')';
                    sens_map{9} = typecast(decode(23:24), 'uint16')';
                    sens_map{10} = typecast(decode(25:26), 'uint16')';

                    % Centre frequency / Sensitivity map pairs
                    % Maximum of 9 pairs frequency / maps stored
                    current_decode_element = 27;
                    for i = 1:sens_map{9}+1 % sens.wnum+1
                        sens_map{11 + 2*(i-1)} = typecast(decode(current_decode_element:current_decode_element+3), 'single')';
                        sens_map{12 + 2*(i-1)} = typecast(decode(current_decode_element+4:(4*double(sens_map{2})*double(sens_map{3})*double(sens_map{4}) + current_decode_element+3)), 'single')'; 
                        current_decode_element = 4*double(sens_map{2})*double(sens_map{3})*double(sens_map{4}) + current_decode_element + 4;
                        % Convert the entries to double to avoid reaching maximum representable number with uint16 data-type
                        % [1 sens_file_content{2}*sens_file_content{3}*sens_file_content{4}] = [1 sens.nx*sens.ny*sens.nz]
                    end

                    % Contour coordinates of the crop map
                    sens_map{13 + 2*(i-1)} = typecast(decode(current_decode_element:current_decode_element + 4*double(sens_map{10})*3 - 1), 'single')';
                    % Convert the entries to double to avoid reaching maximum representable number with uint16 data-type
                    % [1 sens.nb_contour_points*3] = [1 sens_file_content{10}*3]

                    obj.OA_sens_map = sens_map;
                end
            end
                
            % Decoding the US sens map
            if dom.isSetUltrasoundSensitivityMap
                % Because the sensitivity map contains two different data types (both
                % uint16 and single), the base64 decoding has to be performed with care. 
                % Indeed each entry of the file should be decoded according to its own type.
                % The content of the .sens file is stored in a cell array,
                % allowing for mixed data-type (unlike vectors or arrays) with one element per cell.

                code = char(dom.xgetUltrasoundSensitivityMap.getStringValue);
                
                if isempty(code)
                    obj.US_sens_map = {};
                    % Covers the case when the sens map field exists but is empty
                else
                    decode = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(uint8(code)), 'uint8');
                    sens_map = {};

                    % Decode header
                    sens_map{1} = typecast(decode(1:2), 'uint16')';
                    sens_map{2} = typecast(decode(3:4), 'uint16')';
                    sens_map{3} = typecast(decode(5:6), 'uint16')';
                    sens_map{4} = typecast(decode(7:8), 'uint16')';

                    sens_map{5} = typecast(decode(9:12), 'single')';
                    sens_map{6} = typecast(decode(13:16), 'single')';
                    sens_map{7} = typecast(decode(17:20), 'single')';

                    sens_map{8} = typecast(decode(21:22), 'uint16')';
                    sens_map{9} = typecast(decode(23:24), 'uint16')';
                    sens_map{10} = typecast(decode(25:26), 'uint16')';

                    % Centre frequency / Sensitivity map pairs
                    % Maximum of 9 pairs frequency / maps stored
                    current_decode_element = 27;
                    for i = 1:sens_map{9}+1 % sens.wnum+1
                        sens_map{11 + 2*(i-1)} = typecast(decode(current_decode_element:current_decode_element+3), 'single')';
                        sens_map{12 + 2*(i-1)} = typecast(decode(current_decode_element+4:(4*double(sens_map{2})*double(sens_map{3})*double(sens_map{4}) + current_decode_element+3)), 'single')'; 
                        current_decode_element = 4*double(sens_map{2})*double(sens_map{3})*double(sens_map{4}) + current_decode_element + 4;
                        % Convert the entries to double to avoid reaching maximum representable number with uint16 data-type
                        % [1 sens_file_content{2}*sens_file_content{3}*sens_file_content{4}] = [1 sens.nx*sens.ny*sens.nz]
                    end

                    % Contour coordinates of the crop map
                    sens_map{13 + 2*(i-1)} = typecast(decode(current_decode_element:current_decode_element + 4*double(sens_map{10})*3 - 1), 'single')';
                    % Convert the entries to double to avoid reaching maximum representable number with uint16 data-type
                    % [1 sens.nb_contour_points*3] = [1 sens_file_content{10}*3]

                    obj.US_sens_map = sens_map;
                end
            end
        end      
    end
    
    function dom = loadXML(filename)
        % load XMl file in filename to DOM Object
        if ~exist(filename,'file')
            error('File does not exist: %s',filename);
        end
        
        try
            patdoc = javaMethod('parse', 'com.itheramedical.msotbeans.DataModelSensorPositionsDocument$Factory', java.io.File(filename));
            dom = patdoc.getDataModelSensorPositions();
        catch
            try
                patdoc = javaMethod('parse', 'com.itheramedical.msotbeans.DataModelProbeRedistributableDocument$Factory', java.io.File(filename));
                dom = patdoc.getDataModelProbeRedistributable();
            catch ex
                error(['Selected MSOT Probe file could not be loaded. Check well-formedness. ' ex.message]) ;
            end
        end 
    end
    
    function obj = create(di)
        % create Probe Object based on msotData object from Scan
        if ~isvalid(di) || ~isa(di,'msotData'), error('Not a valid msotData object passed to msotProbe.create'); end
        
        % Generic Probe Title
        prid = 'unknownProbe';
        if isempty(di.HWDesc.DeviceSN), did = 'unknown'; else, did = di.HWDesc.DeviceSN; end
        
        % Generate parametric object from old datasets
        if isempty(di.HWDesc.ProjectionData)
           cov = (di.HWDesc.EndAngle - di.HWDesc.StartAngle) / pi * 180;
           obj = msotProbeParametric(prid,did,di.HWDesc.Radius,cov,di.HWDesc.NumDetectors,true);
        % Generate Coordinate based definition
        else
           obj = msotProbeCoordinates(prid,did,di.HWDesc.ProjectionData,true);
        end
        
        % Light Spot Settings
        obj.LightSpotSize = di.HWDesc.LightSpotSize;
        obj.AxialOffset = di.HWDesc.AxialOffset;
        
        % Couplant
        if di.HWDesc.CouplantSOS ~= 0
            obj.CouplantSOS = di.HWDesc.CouplantSOS;
            obj.CouplantCOR = di.HWDesc.CouplantCOR;
            obj.CouplantCurvature = di.HWDesc.CouplantCurvature;
        end
        
        % Impulse Response
        obj.ImpulseResponse = di.loadImpulseResponse;
        
        obj.muted = false;
    end
end

end