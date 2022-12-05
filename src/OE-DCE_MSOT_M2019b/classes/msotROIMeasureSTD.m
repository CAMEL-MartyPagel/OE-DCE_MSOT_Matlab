classdef msotROIMeasureSTD < msotROIMeasure
    properties (Constant)
        Description = 'Standard Deviation (neg=zero)';
        Short = 'STD';
        Ident = 'STD';
        Priority = uint16(1);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val(val < 0) = 0;
            mx = nanstd(val(:));
        end
    end
end
