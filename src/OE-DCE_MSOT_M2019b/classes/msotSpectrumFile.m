classdef msotSpectrumFile < msotSpectrum,

    properties (SetAccess = private),
       wavelengths@uint16 vector = uint16([660:1400]);
       data@single vector = single([]); 

       title@char vector = '';
       filename@char vector = '';
    end

    properties (Constant)
        DIRS = {'C:\ProgramData\iThera\ViewMSOTc\Factory Spectra',...
                'C:\ProgramData\iThera\MSOTSystem\Factory Spectra',...
                'C:\ProgramData\iThera\MSOTSystem\Spectra'};
    end
    
    
    methods
        function obj = msotSpectrumFile(fname,varargin)
            obj@msotSpectrum(varargin);
            
            if exist(fname,'file'), readFile(obj,fname); return; end
            if exist([fname '.csv'],'file'), readFile(obj,[fname '.csv']); return; end
            for j = 1:numel(obj.DIRS),
                if exist([obj.DIRS{j} filesep fname],'file'), readFile(obj,[obj.DIRS{j} filesep fname]); return; end                
                if exist([obj.DIRS{j} filesep fname '.csv'],'file'), readFile(obj,[obj.DIRS{j} filesep fname '.csv']); return; end                
            end
            
            error('Spectrum %s not found',fname);
        end
            
        function readFile(obj,fname)
            obj.filename = fname;
            sfile = importdata(fname,';');

            clear wl abs;
            obj.data = zeros(numel(sfile),1,'single');
            obj.wavelengths = zeros(numel(sfile),1,'uint16');
            for l = 1:numel(sfile)
                tmp = textscan(strrep(sfile{l},',','.'),'%f;%f');
                obj.wavelengths(l) = tmp{1};
                obj.data(l) = tmp{2};
            end
        end
        
    end
end