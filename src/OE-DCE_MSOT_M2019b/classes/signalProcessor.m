classdef signalProcessor < handle
    % signalProcessor processes MSOT Signals in Frequency Domain
    %  Signals are Filtered and/or deconvolved depending on the presence
    %  of msotEIR and msotFilter classes
    properties
        SIG@single = single([]);            % Input Signals (samples x channels) in single
        ImpResp@msotEIR scalar;             % Impulse Repsonse Object (msotEIR)
        Filter@msotFilter vector;           % MSOT Filter Object (msotFilter)
      
        nFFT@uint16 = uint16(4096);         % Length of FFT (default: 4096)
        
        applyFilter@logical = true;         % apply filter in process method
        applyImpResp@logical = true;        % apply ImpulseResponse in process
        progress@logical = false;
    end
        
    methods
        function sig = filter(obj,varargin)
            % filter  Filter the input Signal
            % Filters the input signal using an MSOT Filter. Per default
            % the SIG and Filter properties are used, but may be overridden
            % by the optional parameters
            % 
            % sig = signalProcessor.filter();
            % sig = signalProcessor.filter(sig);
            % sig = signalProcessor.filter(sig, filter);
            sig = getSignal(obj,varargin);
            
            % decide which Filter Object to use
            if nargin >= 3 && isvalid(varargin{2}),
                filt = varargin{2};
            else
                filt = obj.Filter;
            end
            if isempty(filt) || ~isvalid(filt), error('MSOT:signalProcessor:NoFilter','No filter given'); end
            
            sig = filt.filter(sig,obj.progress);
        end
        
        function sig = deconvolve(obj,varargin)
            % deconvolve  Deconvolve the input Signal
            % Deconvolves the input signal using an MSOT Filter. Per default
            % the SIG and ImpResp properties are used, but may be overridden
            % by the optional parameters
            % 
            % sig = signalProcessor.deconvolve();
            % sig = signalProcessor.deconvolve(sig);
            % sig = signalProcessor.deconvolve(sig, eir);
            sig = getSignal(obj,varargin);
            
            % decide which Filter Object to use
            if nargin >= 3 && isvalid(varargin{2}),
                eir = varargin{2};
            else
                eir = obj.ImpResp;
            end
            if isempty(eir) || ~isvalid(eir), error('MSOT:signalProcessor:NoEIR','No impulse Response given'); end
            
            sig = eir.deconvolve(sig,obj.progress);
        end
        
        
        
        
        
        function sig = filter_deconv(obj,varargin)
            % filter_deconvolve  Filter and Deconvolve the input Signal
            % Filters and Deconvolves the input signal using an MSOT Filter. Per default
            % the SIG, Filter and ImpResp properties are used, but may be overridden
            % by the optional parameters
            % 
            % sig = signalProcessor.deconvolve();
            % sig = signalProcessor.deconvolve(sig);
            % sig = signalProcessor.deconvolve(sig, filt);
            % sig = signalProcessor.deconvolve(sig, filt, eir);
            sig = getSignal(obj,varargin);
            svec = size(sig);
            nsamples = svec(1);
            
            wbar = [];
            if obj.progress, wbar = waitbar(0,'Filtering and Deconvolution...'); end;

            % transform to frequency domain if necessary
            if isreal(sig),
                if obj.progress, waitbar(0,wbar,'FFT...'); end;
                sig = sig - repmat(mean(sig,1),[nsamples 1]);
                sig = fft(sig,obj.nFFT,1);
                tdom = true;
            else
                tdom = false;
            end
            
            % decide which EIR Object to use
            if nargin >= 4 && isvalid(varargin{3}),
                eir = varargin{3};
            else
                eir = obj.ImpResp;
            end
            
            % Deconvolve, if appropriate
            if ~isempty(eir) && isvalid(eir), 
                sig = eir.deconvolve(sig,wbar);
            end
            
            
            % decide which Filter Object to use
            if nargin >= 3 && isvalid(varargin{2}),
                filt = varargin{2};
            else
                filt = obj.Filter;
            end
            
            % Filter, if appropriate
            if ~isempty(filt) 
                for jj = 1:numel(filt),
                    if isvalid(filt(jj)), sig = filt(jj).filter(sig,wbar); end
                end
            end

            % transform back to time domain
            if tdom,
                if obj.progress, waitbar(0,wbar,'iFFT...'); end;
                sig = real(ifft(sig,[],1));
                sig = reshape(sig(1:nsamples,:),svec);
            end
            
            if obj.progress, close(wbar); end;
            
        end
        
        
        
        
        function sig = process(obj,sig)
            % process  Apply selected processing methods
            filt = false;
            if ~isempty(obj.Filter) && isvalid(obj.Filter(1)) && obj.applyFilter, filt = true; end
            
            eir = false;
            if ~isempty(obj.ImpResp) && isvalid(obj.ImpResp) && obj.applyImpResp, eir = true; end
            
            if filt && eir,
                sig = filter_deconv(obj,sig);
            elseif filt
                sig = filter(obj,sig);
            elseif eir
                sig = deconvolve(obj,sig);
            end
        end
        
    end
    
    
    
    methods (Access = private)
        
        % *** Auxiliaries
        function sig = getSignal(obj,args)
            if nargin >= 2 && ~isempty(args{1}),
                sig = args{1};
            else
                sig = obj.SIG;
            end
            if isempty(sig), error('MSOT:signalProcessor:NoInputSignal','No input signal given'); end
        end
        
        
    end
end