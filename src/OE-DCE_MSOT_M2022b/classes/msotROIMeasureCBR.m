classdef msotROIMeasureCBR < msotROIMeasure
    properties (Constant)
        Description = 'Contrast-to-Background Ratio';
        Short = 'CBR';
        Ident = 'CBR';
        Priority = uint16(0);
    end
    
    methods
        function mx = calculate(obj,~,img,roi)
            % Foreground Values
            val = img(roi);
            val(val < 0) = 0;
            
            % Background Pixels
            valbg = img(~roi);
            
            % CBR (ignoring NaNs - these were likely intentionally
            % filtered)
            mx = nanmean(val) ./ nanstd(valbg);
        end
    end
end
