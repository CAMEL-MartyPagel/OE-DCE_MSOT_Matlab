classdef msotProbeParametric < msotProbe
    % msotProbeParametric Parametric MSOT Probe Definition
    %   Parametric Definition of an MSOT Detector, determined by Radius,
    %   Coverage and N
    % 
    % see help for msotProbe
    
properties (Constant = true)
    Type@char vector = 'Parametric';
    GeometryDefinitionType@char vector = 'Parametric';
    % Two properties coexist to maintain field names of the legacy scheme
end

properties
    Radius@double = 0;    % radius of detector (only Parametric)
    Coverage@double = 0;
    NumDetectors@uint16 = uint16(0);
end

properties (Dependent = true)
    Angles@double;            % angle vector, will determine Sensors, N and Coverage when set directly
    StartAngle@double;        % Angle of first Element, read only
    StepAngle@double;         % Step Angle, read only 
    EndAngle@double;          % End Angle, read only 
    Pitch@double;             % Width of an Element, in Plane (only Parametric)
end

properties (Dependent = true, SetObservable = true)
    Sensors@double;           % Sensors vector, can only be set for Type = Coordinates, otherwise is automatically calculated
end    
  
methods
    % *****************************************************************
    %% CONSTRUCTOR
    function obj = msotProbeParametric(ID, DID, R, C, Ne, varargin)
        % Create MSOT Parametric Probe Object
        %
        % msotProbeParametric(ProbeID,DesignID,Radius,Coverage,N)
        obj@msotProbe(ID,DID);
        obj.Radius = R;
        obj.Coverage = C;
        obj.NumDetectors = uint16(Ne);
        
        % if coordinates given (if read from File) validate Geometry
        if numel(varargin) >= 1
            S_ = varargin{1};
            
            obj.Handheld = varargin{2};
            
            if numel(varargin) >= 3
                obj.ViewMSOTv4Compatible = varargin{3};
            end
            
            S = obj.Sensors;
            if numel(S) ~= numel(S_)
                error('Inconsistent Probe File - Number of Sensors does not match! (%i - %i)',numel(S)/3,numel(S_)/3);
            end
            Sdiff = abs(S-S_);
            if max(Sdiff(:)) > eps
                for jj = 1:size(S,1)
                    fprintf('%.17f\t%.17f\t%.17f\n',Sdiff(jj,:));
                end
                warning('Inconsistent Probe File - Coordinates dont match description!');
            end
        end
    end
    
    
    % *****************************************************************
    %% ACCESSORS
    function StepAngle = get.StepAngle(obj)
        StepAngle = ( obj.Coverage * (pi/180) ) / double(obj.N); 
    end

    function StartAngle = get.StartAngle(obj)
        StartAngle = obj.Angles(1);
    end

    function EndAngle = get.EndAngle(obj)
        EndAngle = obj.Angles(end);
    end
    
    function A = get.Angles(obj) 
        SA = -0.5*obj.Coverage*(pi/180) + obj.StepAngle/2  - pi/2; 
        if ~obj.Handheld 
            SA = SA + pi; 
        end
        A = SA+double((1:obj.N/2)-1)*obj.StepAngle + pi/2;
        A = unwrap([A -fliplr(A)]-pi/2);
        if obj.ViewMSOTv4Compatible
            A = A - pi;
        end
    end
    
    function pos_elements = get.Sensors(obj)
        pos_elements=zeros(obj.N,3);
        if ~obj.ViewMSOTv4Compatible
            [pos_elements(:,1),pos_elements(:,3),pos_elements(:,2)]=sph2cart(zeros(obj.N,1),obj.Angles',ones(obj.N,1)*obj.Radius);
        else  % in the new case, probe orientation is different, where the depth axis is along z
            [pos_elements(:,1),pos_elements(:,2),pos_elements(:,3)]=sph2cart(zeros(obj.N,1),obj.Angles',ones(obj.N,1)*obj.Radius);
        end
    end
    
    function p = get.Pitch(obj)
        S = obj.Sensors;
        p = sqrt((S(1,1)-S(2,1)).^2+(S(1,2)-S(2,2)).^2);
    end
    
    % *****************************************************************
    %% FUNCTIONALITY
    function str = writeXML_(obj)
        str = {};
        str{end+1} = sprintf('\t<Radius>%.17f</Radius>\n',obj.Radius);
        str{end+1} = sprintf('\t<Coverage>%.17f</Coverage>\n',obj.Coverage);
        str = [str{:}];
    end   
    
    function obj = readXML_(obj,dom)
        if dom.isSetNumDetectors
            obj.NumDetectors = dom.getNumDetectors; 
        else
            error('Cannot read Parametric Definition without NumDetectors being set');
        end

        if dom.isSetRadius
            obj.Radius = double(dom.getRadius); 
        else
            error('Cannot read Parametric Definition without Radius being set');
        end
    
        if dom.isSetCoverage
            obj.Coverage = double(dom.getCoverage); 
        else
            error('Cannot read Parametric Definition without Coverage being set');
        end
    end
        
end
end