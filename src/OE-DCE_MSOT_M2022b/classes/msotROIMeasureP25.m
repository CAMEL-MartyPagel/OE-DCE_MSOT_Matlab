classdef msotROIMeasureP25 < msotROIMeasure
    properties (Constant)
        Description = 'Percentile 25';
        Short = 'P25';
        Ident = 'P25';
        Priority = uint16(1);
        Percentile = 25;
    end
    
    methods
        function mx = calculate(obj,val,varargin)
            val(val < 0) = 0;
            val = sort(val(:));
            idx = round(numel(val)*obj.Percentile/100);
            if (idx == 0),idx = 1; end
            mx = val(idx);
        end
    end
end
