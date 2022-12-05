classdef msotROIMeasureMed < msotROIMeasure
    properties (Constant)
        Description = 'Median';
        Short = 'Med';
        Ident = 'Med';
        Priority = uint16(2);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val(val < 0) = 0;
            mx = nanmedian(val(:));
        end
    end
end
