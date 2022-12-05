classdef msotSettingsChangedEventData < event.EventData
    % msotSettingsChangedEventData  Event Object to carry information concerning Settings Modifications
    
    properties
        sID@char vector = '';               % sID of Settings Object
        oldHash@char vector = '';           % hashCode before Modification
        newHash@char vector = '';           % hashCode after Modification
        fieldName@char vector = '';         % Name of changed property
        Description@char vector = '';       % Descrition of Modification (optional)
    end
    
    methods
        function obj = msotSettingsChangedEventData(sid,oH,nH,varargin)
            obj.sID = sid;
            obj.oldHash = oH;
            obj.newHash = nH;
            if numel(varargin) >= 1
                obj.fieldName = varargin{1};
            end
            if numel(varargin) >= 2
                obj.Description = varargin{2};
            end
        end
    end
end