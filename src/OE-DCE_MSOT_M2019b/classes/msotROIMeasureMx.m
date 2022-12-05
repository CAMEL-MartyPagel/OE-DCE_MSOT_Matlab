classdef msotROIMeasureMx < msotROIMeasure
    properties (Constant)
        Description = 'Maximum';
        Short = 'Max';
        Ident = 'Mx';
        Priority = uint16(3);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val(val < 0) = 0;
            mx = nanmax(val(:));
        end
    end
end
