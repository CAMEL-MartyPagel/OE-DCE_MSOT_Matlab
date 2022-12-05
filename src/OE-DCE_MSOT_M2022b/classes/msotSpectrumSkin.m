classdef msotSpectrumSkin < msotSpectrum,

% Based on data from:
% Optical properties of human skin, subcutaneous and mucous tissues in the
% wavelength range from 400 to 2000nm 
% A N Bashkatov1, E A Genina, V I Kochubey and V V Tuchin
% J. Phys. D: Appl. Phys. 38 (2005) 2543–2555

% UNIT: 1/cm


    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([ 
598.88
650.42
700.38
752.04
802.06
847.23
898.83
950.49
998.87
1047.20
1100.46
1150.81
1199.40
1247.62
1302.59
1349.83
1400.84
]);

data@single vector = single([ 
0.69
0.57
0.48
0.47
0.44
0.40
0.34
0.33
0.28
0.19
0.16
0.42
0.56
0.36
0.42
0.79
1.65

]); 

       title@char vector = 'Skin';
    end
   
    methods
        function obj = msotSpectrumSkin(varargin),
            obj@msotSpectrum(varargin);
        end
    end
end