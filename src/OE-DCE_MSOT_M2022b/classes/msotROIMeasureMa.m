classdef msotROIMeasureMa < msotROIMeasure
    properties (Constant)
        Description = 'Mean (all pixels)';
        Short = 'Mall';
        Ident = 'Ma';
        Priority = uint16(7);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val = val(:);
            mx = mean(val);
        end
    end
end
