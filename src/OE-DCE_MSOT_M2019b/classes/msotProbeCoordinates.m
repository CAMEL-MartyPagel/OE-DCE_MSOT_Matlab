classdef msotProbeCoordinates < msotProbe
    
    properties (Constant = true)
        Type@char vector = 'Coordinates';
        GeometryDefinitionType@char vector = 'Coordinates';
        % Two properties coexist to maintain field names of the legacy scheme
    end

    properties (SetObservable)
        Sensors@double = zeros(0,3);
    end
    
    properties
        StartAngle@double = nan;                % Angle of first Element, read only
        StepAngle@double = nan;                 % Step Angle, read only 
        EndAngle@double = nan;                  % End Angle, read only 
        Pitch@double = nan;                     % Width of an Element, in Plane (only Parametric)
    end
    
    properties (Dependent = true)
        Radius@double;
        Coverage@double;
        NumDetectors@uint16;
        Angles@double;
        Spherical@double;
    end
  
methods
    % *****************************************************************
    %% CONSTRUCTOR
    function obj = msotProbeCoordinates(ID, DID, S, varargin)
        % Create MSOT Coordinate Probe Object
        %
        % msotProbeParametric(ProbeID,DesignID,Radius,Coverage,N)
        obj@msotProbe(ID,DID);
        obj.Sensors = S;
        
        if numel(varargin) >= 1
            obj.ViewMSOTv4Compatible = logical(varargin{1});
        end
    end
    
    
    % *****************************************************************
    %% ACCESSORS
    function R = get.Radius(obj)
        if ~obj.muted 
            warning('Radius may be inaccurate for Probes with Coordinate Definition'); 
        end
        R = mean(obj.Spherical(:,3));
    end

    function C = get.Coverage(obj)
        if isempty(obj.Spherical) 
            C = 0; 
            return 
        end
        
        if ~obj.muted 
            warning('Coverage may be inaccurate for Probes with Coordinate Definition'); 
        end
        
        ang_=obj.Spherical(:,2);
        
        % Definition for the Coverage --> Consider the whole cup including the hole on top
        C = round(2 * abs(diff([min(ang_),pi/2])) / pi * 180);
        
        % Alternative definition: Consider only the angle covered by the detectors (limited on top by the hole)
        % C = round(abs(diff([min(ang_),max(ang_)])) / pi * 180);
                
                % Old definition
                % C = C*(obj.N+1)/obj.N; % Gives obviously wrong results
    end

    function N = get.NumDetectors(obj)
        N = uint16(size(obj.Sensors,1));
    end

    function A = get.Angles(obj) 
        if ~obj.muted 
            warning('Angles may be inaccurate for Probes with Coordinate Definition') 
        end        
        A = obj.Spherical(:,2);
    end
       
    function S = get.Spherical(obj)
        Sens = obj.Sensors;
        if obj.ViewMSOTv4Compatible
            [S(:,1), S(:,2), S(:,3)] = cart2sph(Sens(:,1),Sens(:,2),Sens(:,3)); % phi - theta - R
        else
            [S(:,1), S(:,2), S(:,3)] = cart2sph(Sens(:,1),Sens(:,3),Sens(:,2)); % phi - theta - R (adapted to the y-z swap to output always phi-theta-R)
        end

        loc = S(:,1) == pi;
        S(loc,1) = 0;
        % S(loc,2) = -(pi+S(loc,2)); % Gives obviously wrong results
    end
    
    % Overwriting the set method for Sensors property
    function set.Sensors(obj,Sensors)
        if isempty(Sensors) 
            obj.Sensors = []; 
            return 
        end
        if size(Sensors,2) ~= 3
            error('Sensor Position Array does not have three dimensions');
        end

        obj.Sensors = Sensors;
        if ~obj.isSymmetric % Dot access similar here to the more logical isSymmetric(obj), method of the msotProbe class!
            warning('Provided Geometry is not symmetric');
        end
    end
    
    
    % *****************************************************************
    %% FUNCTIONALITY
    function str = writeXML_(obj)
        str = {};
        str{end+1} = sprintf('\t<Radius>%.17f</Radius>\n',obj.Radius);
        str{end+1} = sprintf('\t<Coverage>%.17f</Coverage>\n',obj.Coverage);
        str = [str{:}];
        
        % str = [];
    end
    
    function obj = readXML_(obj,~)
    end
  end
end