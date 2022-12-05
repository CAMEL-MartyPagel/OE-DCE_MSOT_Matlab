classdef msotFilter < handle
    % msotFilter  Abstract Filter Class for MSOT Filters
    properties
        FS@double = 4e7;                  % Sampling Frequency
        nFFT@uint16 = uint16(4096);       % Length of FFT
    end
    
    methods
        
        function fd = fdom(~)
            % Frequency Domain Response (Length defined by nFFT)
            fd = [];
        end
        
        function str = toString(~)
            % string Representation
            str = 'custom';
        end
        
        function configure(~)
        end
 
        function sig_ = filter(obj,sig_,varargin)
            % filter  Filter the input Signal
            % Apply the filter to the given input signal. Signal can be in
            % Time Domain or in Frequency Domain
            % sig_f = filter(sig);
            %
            % Parameters:
            % 1: Signal (samples x channels) 
            if isempty(sig_), error('MSOT:signalProcessor:NoInputSignal','No input signal given'); end
            
            % waitbar
            wbar = [];
            if numel(varargin) >= 1,
                if isa(varargin{1},'logical') && varargin{1}, wbar = waitbar(0,'Filtering...'); 
                elseif isa(varargin{1},'matlab.ui.Figure'), wbar = varargin{1}; end
            end

            svec = size(sig_);
            nsamples = svec(1);
            nchannels = svec(2);

            % transform to frequency domain if necessary
            if isreal(sig_),
                if ~isempty(wbar), waitbar(0,wbar,'FFT...'); end;
                sig_ = sig_ - repmat(mean(sig_,1),nsamples,1);
                sig_ = fft(sig_,obj.nFFT,1);
                tdom = true;
            else
                tdom = false;
            end

            % filter response
            GG = repmat(single(obj.fdom), 1, nchannels);

            sig_ = reshape(sig_,[size(GG) prod(svec(3:end))]);
            for jj = 1:size(sig_,3),
                % apply filter
                if ~isempty(wbar), waitbar(jj/size(sig_,3),wbar,'Filtering...'); end;
                sig_(:,:,jj) = sig_(:,:,jj) .* GG;
            end
                            % transform back to time domain
            if tdom,
                if ~isempty(wbar), waitbar(0,wbar,'iFFT...'); end;
                sig_ = real(ifft(sig_,[],1));
                sig_ = sig_(1:nsamples,:);
            end

            sig_ = reshape(sig_,svec);
            if numel(varargin) >= 1 && isa(varargin{1},'logical') && varargin{1}, close(wbar); end
        end
    end        
end