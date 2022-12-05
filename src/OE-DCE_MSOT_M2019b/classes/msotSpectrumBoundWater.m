classdef msotSpectrumBoundWater < msotSpectrum,
% Chung, So Hyun, Hon Yu, Min-Ying Su, Albert E. Cerussi, and Bruce J. Tromberg. ‘Molecular Imaging of Water Binding State and Diffusion in Breast Cancer Using Diffuse Optical Spectroscopy and Diffusion Weighted MRI’. Journal of Biomedical Optics 17, no. 7 (2012): 071304–071304. doi:10.1117/1.JBO.17.7.071304.
    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([ 900
910
920
930
940
950
960
970
980
990
1000]);

data@single vector = single([ 0
0
0
0
0.08
0.2
0.3
0.41
0.42
0.41
0.4]);

       title@char vector = 'BoundWater';
   end
   
   methods
        function obj = msotSpectrumBoundWater(varargin)
            obj@msotSpectrum(varargin);
        end
   end

end