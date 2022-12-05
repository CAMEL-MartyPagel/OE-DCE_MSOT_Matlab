classdef msotROIMeasurePPP < msotROIMeasure
    properties (Constant)
        Description = 'Percent Positive Pixels';
        Short = 'PPP';
        Ident = 'PPP';
        Priority = uint16(0);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            mx = nnz(val(:) >= 0) / numel(val) * 100;
        end
    end
end
