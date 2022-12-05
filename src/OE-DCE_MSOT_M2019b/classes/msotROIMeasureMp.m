classdef msotROIMeasureMp < msotROIMeasure
    properties (Constant)
        Description = 'Mean (only positives)';
        Short = 'Mpos';
        Ident = 'Mp';
        Priority = uint16(9);
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val = val(:);
            val = val(val >= 0);
            mx = nanmean(val);
        end
    end
end
