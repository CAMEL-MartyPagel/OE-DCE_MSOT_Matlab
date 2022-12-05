classdef msotLaserCalibrationElement < handle
    properties
        Wavelength@uint16 scalar = uint16(0);
        PulseDivider@uint16 scalar = uint16(1);
        DividerValue@uint16 scalar = uint16(nan);       % Alias for Compatibility
        AverageEnergy@double scalar = 0;
        Attenuation@double scalar = 0;
        Coefficients@double vector = [];
        Coefficiants@double vector = [];                % Alias for Compatibility
    end
    
    methods
        function obj = msotLaserCalibrationElement(wl,pd,ae,att,co)
            if nargin == 0, return; end
            
            obj.Wavelength = wl; 
            obj.PulseDivider = pd; 
            obj.AverageEnergy = ae; 
            obj.Attenuation = att;
            if nargin >= 5, obj.Coefficients = co; end
        end
        
        function en = calculate(obj,diode)
            % Calculate Energy Value based on supplied DiodeReadout
            en = polyval(flip(obj.Coefficients),diode);
        end
        
        function ret = get.Coefficiants(obj)
            % Alias for misspelled Attribute in original code to maintain
            % compatibility
            ret = obj.Coefficients;
        end
        
        function ret = get.DividerValue(obj)
            % Alias for Attribute in original code to maintain compatibility
            ret = obj.PulseDivider;
        end
        
    end
end