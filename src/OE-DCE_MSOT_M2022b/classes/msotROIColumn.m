classdef msotROIColumn
    properties
        layer@char vector = char([]);
        measure = msotROIMeasureMz();
    end
    
    properties (Constant),
        measureSep = '|';
        arraySep = ',';
    end
    
    methods
        function obj = msotROIColumn(lay,meas)
            if nargin == 0, return;
            elseif nargin == 1,
                [obj.layer, obj.measure] = fromString(obj,lay);
            elseif nargin == 2,
                obj.layer = lay;
                obj.measure = meas;
            end
        end
        
        
        % *********************************************************************
        % Accessor functions
        function obj = set.measure(obj,str)
            if (isa(str,'msotROIMeasure')); obj.measure = str; 
            elseif ischar(str),
                obj.measure = msotROIMeasure.create(str);
            else
                error('Value for Property measure must be char or msotROIMeasure');
            end
        end
        
        
        
        % *********************************************************************
        % String conversion
        function [lay, meas] = fromString(obj,str)
            if iscell(str), str = char(str); end
            tok = strsplit(str,obj.measureSep);
            if isempty(tok), 
                error('Cannot parse %s',str);
            end
            if numel(tok) == 1,
                error('Missing input for measure (Format: "layer%smeasure")',obj.measureSep);
            end
            lay = tok{1};
            meas = tok{2};
        end
        
        function str = toString(obj)
            str = [obj.layer obj.measureSep obj.measure.Ident];
        end
    end
    
    
    
    % *********************************************************************
    % Static Properties to create Arrays
    methods (Static)
        function arr = parseString(str)
            str = strsplit(str,msotROIColumn.arraySep);
            arr = msotROIColumn.empty(0,numel(str));
            for jj = 1:numel(str)
                arr(jj) = msotROIColumn(str{jj});
            end
        end
        
        function str = createString(arr)
            if ~isa(arr,'msotROIColumn'), error('Wrong input type (expecting msotROIColumn)'); end
            slist = cell(1,numel(arr));
            for jj = 1:numel(arr),
                slist{jj} = arr(jj).toString();
            end
            str = strjoin(slist,msotROIColumn.arraySep);
        end
    end
end