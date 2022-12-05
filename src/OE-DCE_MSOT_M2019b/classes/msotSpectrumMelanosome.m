classdef msotSpectrumMelanosome < msotSpectrum,
% based on melanosomes
%
% References: 
% -----------
% L Goldman, The Skin, Arch. Environmental Health, 18:435, 1969.
% DH Sliney, WA Palmisano. The evaluation of laser hazards. AIHA Journ. 20:425, 1968.
% SL Jacques, DJ McAuliffe. The melanosome: threshold temperature for explosive vaporization and internal absorption coefficient during pulsed laser irradiation. Photochem. Photobiol. 53:769-775, 1991.
% SL Jacques, RD Glickman, JA Schwartz. Internal absorption coefficient and threshold for pulsed laser disruption of melanosomes isolated from retinal pigment epithelium. SPIE Proceedings 2681:468-477, 1996.
 
    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([660:1400])';
       data@single vector = single([]); 

       title@char vector = 'Melanin';
    end
   
    methods
        function obj = msotSpectrumMelanosome(varargin),
            obj@msotSpectrum(varargin);
        end
        
        function val = get.data(obj)
            val = single(1.7e12*(double(obj.wavelengths).^-3.48));
        end
    end
end