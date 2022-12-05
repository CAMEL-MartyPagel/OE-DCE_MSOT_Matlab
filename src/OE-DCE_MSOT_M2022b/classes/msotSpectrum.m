classdef (Abstract) msotSpectrum < handle
    % msotSpectrum  Spectra base class
    %
    % To form multiples use a cell array, and use 
    %    cell2mat(cellfun(@(x) x(700:10:900),arr','UniformOutput',false))
    % to extract spectra as matrix for specified wavelengths
    properties (Abstract,SetAccess = private)
        data@single vector;
        wavelengths@uint16 vector;
        title@char vector;
    end
    
    properties
        WL@uint16;
    end
    
    methods
        function obj = msotSpectrum(wls)
            if nargin == 0, return; end
            if iscell(wls), 
                if ~isempty(wls), wls = wls{1}; 
                else wls = []; end;
            end
            obj.WL = uint16(wls);
        end
        
        function B = subsref(obj,S)
            A = S.subs{1};
            % convert wavelengths to indices
            if min(A) > 300,
                [locA,locB] = ismember(obj.wavelengths,A);
                B = zeros(size(A),'single');
                DONE = []; locB = locB(locB > 0);
                % matching items
                if nnz(locA)>0, 
                    B(locB) = obj.data(locA); 
                    DONE = A(locB);     % mark those as done
                end
                [locDONE] = ismember(A,DONE); % logical version of B...
                locI = ~locDONE & (A > min(obj.wavelengths)) & (A < max(obj.wavelengths));  
                if nnz(locI) > 0,
                    B(locI) = interp1(single(obj.wavelengths),obj.data,A(locI),'pchip',0);
                end
                % rest remains zero
            % use indices
            else
                if isempty(obj.WL),
                    B = obj.data(A);
                else
                    locA = ismember(obj.wavelengths,obj.WL(A));
                    if nnz(locA) == numel(obj.WL(A)),
                        B = obj.data(locA);
                    else
                        if strcmp(A,':'),
                            A = obj.WL;
                        end
                        B = interp1(single(obj.wavelengths),obj.data,single(A),'pchip',0)';
                    end
                end
            end                
           
        end
    end
    
    methods (Static)
        function arr = create(specs,wls)
            if nargin == 1, wls = []; end
            for jj = 1:numel(specs),
                try
                    arr{jj} = eval(['msotSpectrum' specs{jj} '(wls)']);
                catch ex
                    if strcmp(ex.identifier,'MATLAB:UndefinedFunction')
                        arr{jj} = msotSpectrumFile(specs{jj},wls);
                    else
                        rethrow(ex);
                    end
                end
            end
        end
    end
    
end
            