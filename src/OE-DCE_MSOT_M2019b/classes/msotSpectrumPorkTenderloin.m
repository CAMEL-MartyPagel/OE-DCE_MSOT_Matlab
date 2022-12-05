classdef msotSpectrumPorkTenderloin < msotSpectrum,

% Based on data from:
% Optical absorption and scattering properties of bulk porcine muscle
% phantoms from interstitial radiance measurements in 650–900 nm range
% Serge Grabtchak, Logan G Montgomery and William M Whelan
% Phys. Med. Biol. 59 (2014) 2431–2444

% UNIT: 1/cm


    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([ 
610.27
611.60
613.92
621.54
633.81
649.05
668.61
685.84
707.72
726.28
744.51
759.75
776.99
796.21
816.10
840.96
860.84
882.05
900.28
922.16
936.08
947.02
]);

data@single vector = single([ 
0.60
0.56
0.48
0.37
0.29
0.24
0.20
0.17
0.15
0.14
0.14
0.15
0.14
0.12
0.12
0.12
0.12
0.13
0.14
0.15
0.18
0.21

]); 

       title@char vector = 'Skin';
    end
   
    methods
        function obj = msotSpectrumPorkTenderloin(varargin),
            obj@msotSpectrum(varargin);
        end
    end
end