classdef msotSpectrumSubcutis < msotSpectrum,

% Based on data from:
% Optical properties of human skin, subcutaneous and mucous tissues in the
% wavelength range from 400 to 2000nm 
% A N Bashkatov1, E A Genina, V I Kochubey and V V Tuchin
% J. Phys. D: Appl. Phys. 38 (2005) 2543–2555

% UNIT: 1/cm


    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([ 
566.51
620.78
676.88
728.28
798.42
843.63
899.78
951.25
1002.67
1049.44
1111.75
1171.05
1225.57
1278.45
1333.02
1390.88
1445.66
1498.5
]);

data@single vector = single([ 
1.39
1.16
1.13
1.09
1.07
1.06
1.07
1.08
1.06
1.05
1
1.04
1
0.9
0.89
1.02
1.18
1.05
]); 

       title@char vector = 'Skin';
    end
   
    methods
        function obj = msotSpectrumSubcutis(varargin),
            obj@msotSpectrum(varargin);
        end
    end
end