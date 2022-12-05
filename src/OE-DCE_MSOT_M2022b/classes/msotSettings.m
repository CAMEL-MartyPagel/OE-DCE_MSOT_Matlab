classdef (Abstract) msotSettings < msotHandle
    % msotSettings Abstract base class for msot Settings objects
    %  Store Settings, including the ability to generate a uniquely
    %  identifying string, conversion to json and saving
    properties (Constant, Abstract)
        version@uint8;
        prefix@char scalar;
    end
    
    properties (Constant)
        % properties that are not taken into account for hash calculation
        % and detection of changes
        ignoredProperties = {'version','hashCode','ID','sID','logger','creation','undoProp','undoValue','muted'};
    end
    
    properties (SetAccess = protected)
        hashCode@char vector = '';
        undoing@logical = false;
        undoProp@char vector = '';
        undoValue;
    end
    
    properties
        muted@logical = false;      % Mute notification of property changes
    end
    
    events
        Modified;
    end
    
    methods
        % CONSTRUCTOR
        function obj = msotSettings(varargin)
            obj@msotHandle();   % superclass constructor
            
            % set muting if provided
            if numel(varargin) >= 1 && ~isempty(varargin{1})
                obj.muted = logical(varargin{1});
            end
            
            % Set initial Hash Code (only defaults)
            obj.hashCode = generateHashCode(obj);
            obj.logger.debug(obj,'Initial Hash Code: %s',obj.hashCode);
            
            % Associate Listeners to track property changes
            createListeners(obj);
        end
        
        function undo(obj)
            % undo last set action (e.g. if value is invalid)
            if isempty(obj.undoProp) 
                obj.logger.error(obj,'msotSettings:noUndoPossible','No action recorded that can be undone');
                return;
            end
            obj.undoing = true;
            obj.(obj.undoProp) = obj.undoValue;
            obj.hashCode = obj.generateHashCode();
            obj.logger.info(obj,'Undone Change to Property %s - hash now %s',obj.undoProp, obj.hashCode);
            obj.undoProp = '';
            obj.undoing = false;
        end
        
        function str = char(obj,varargin)
            % Generate Property String; Converter to char
            
            % Struct to encode
            prop = struct();
            mc = metaclass(obj);
            % if called as an actual converter, also add version and hash
            % and non-observable properties for tracebility. 
            % Otherwise don't include (for hash calculation)
            includeNonObservable = false;
            if ~(numel(varargin) >= 1 && varargin{1})
                includeNonObservable = true;
                prop.class = mc.Name;
                prop.version = obj.version;
                prop.hashCode = obj.hashCode;
            end
            
            proplist = mc.PropertyList;
            ind = find(~ismember({proplist.Name},[obj.ignoredProperties 'ignoredProperties']));
            for ix = ind
                pname = proplist(ix).Name;
                
                % put this in try/catch just to better find out about weird
                % issues in conversion to char
                try
                    % only observable properties are relevant (they should
                    % influence the hash)
                    if ~includeNonObservable && ~proplist(ix).SetObservable, continue; end
                    % ignore dependent properties only if not obervable
                    if proplist(ix).Dependent && ~proplist(ix).SetObservable, continue; end
                    % always ignore constant or transient properties
                    if proplist(ix).Constant || proplist(ix).Transient, continue; end 
                    
                    pval = obj.(pname);
                    % ignore properties with empty default and value
                    if proplist(ix).HasDefault && isempty(proplist(ix).DefaultValue) && isempty(pval), continue; end
                    % ignore properties set to their default
                    if proplist(ix).HasDefault && numel(proplist(ix).DefaultValue) == numel(pval) && all(all(proplist(ix).DefaultValue == pval)), continue; end
                
                catch ex
                    obj.logger.error(obj,'msotSettings:determineCharError','Error determining whether to include property %s', pname, ex);
                end
                % if the property actually is another class try to convert
                % to char
                if ~isempty(pval) && numel(pval) == 1 && ishandle(pval)
                    try 
                        pval = char(pval);
                    catch err
                        obj.logger.warn(obj,'msotSettings:cannotConvertToChar','Cannot convert property %s to char (%s) - skipping',pname,err.message);
                        continue;
                    end
                end
                        
                % add to property struct
                prop.(pname) = pval;
            end
%             str = jsonencode(prop);       % MATLAB integrated, saves un-pretty output
            str = savejson('',prop,struct('ParseLogical',1));
        end
        
        % override parent to include prefix
        function hash = generateHashCode(obj)
            hash = [obj.prefix '-' generateHashCode@msotHandle(obj)];
        end
              
        function set.muted(obj,m)
            if ~m
                obj.hashCode = obj.generateHashCode;
                obj.logger.info(obj,'Object unmuted, Modifications are now tracked (current hash: %s)',obj.hashCode);
            end
            obj.muted = m;
        end
            
    end
    
    methods (Access = protected)
        function createListeners(obj)
            % Associate all events with PostSet listeners to re-calculate
            % hash and notify potential Listeners
            mc = metaclass(obj);
            proplist = mc.PropertyList;
            ind = find(~ismember({proplist.Name},[obj.ignoredProperties 'ignoredProperties']));
            for ix = ind
                % Ignore constants and not observable properties (will fail...)
                if proplist(ix).Constant || ~proplist(ix).SetObservable, continue; end 
                addlistener(obj,proplist(ix).Name,'PreSet',@obj.propertyNotYetModified);
                addlistener(obj,proplist(ix).Name,'PostSet',@obj.propertyModified);
            end
        end
    end
    
    methods (Access = protected, Static)
        function propertyModified(src,evt)
            % ignore if muted
            if evt.AffectedObject.muted, return; end
            % Trigger the Modified event as soon as a property is changed
            hC = generateHashCode(evt.AffectedObject);
            % If a change actually happened
            if ~strcmp(hC,evt.AffectedObject.hashCode)
                msg = sprintf('Property %s modified - hash now %s',src.Name,hC);
                evt.AffectedObject.logger.debug(evt.AffectedObject,msg);      % debug just in case (with old hash code)
                newevt = msotSettingsChangedEventData(evt.AffectedObject.sID,evt.AffectedObject.hashCode,hC,src.Name,msg);
                evt.AffectedObject.hashCode = hC;
                notify(evt.AffectedObject,'Modified',newevt);
            end
        end
            
        function propertyNotYetModified(src,evt)
            % listener callback just befor setting property
            
            % Only save value if not in course of an undo action and if
            % there are listeners to the Modified object (we want to avoid
            % accidentally undoing changes from a msotProcessor that have been
            % done before creation of the msotProcessor)
            if ~evt.AffectedObject.muted && ~evt.AffectedObject.undoing && event.hasListener(evt.AffectedObject,'Modified')
                evt.AffectedObject.logger.debug(evt.AffectedObject,'Property %s about to be modified - saving value',src.Name);
                evt.AffectedObject.undoProp = src.Name;
                evt.AffectedObject.undoValue = evt.AffectedObject.(src.Name);
            end
        end
    end
end