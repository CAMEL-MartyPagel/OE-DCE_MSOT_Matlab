classdef msotROIMeasureMz < msotROIMeasure
    properties (Constant)
        Description = 'Mean (count neg as zero)';
        Short = 'Mzero';
        Ident = 'Mz';
        Priority = uint16(10);
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val = val(:);
            val(val < 0) = 0;
            mx = nanmean(val);
        end
    end
end
