classdef msotROIMeasureMi < msotROIMeasure
    properties (Constant)
        Description = 'Minimum';
        Short = 'Min';
        Ident = 'Mi';
        Priority = uint16(4);        
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val(val < 0) = 0;
            mx = nanmin(val(:));
        end
    end
end
