classdef msotEIR < handle
    % msotEIR Electrical Impulse Response of MSOT Transducer
    %   Can be used to generate a custom electrical impulse response of a
    %   transducer of the given properties. The impulse response can be
    %   retrieved in frequency and time domain, and also saved to file
    
    properties
      f_center@double = 5e6;            % Center Frequency in Hz (default: 5MHz)
      bw@double = 0.6;                  % Relative T/R Bandwidth (default: 0.6)
      phase_shift@double = 0.0;         % Phase Shift (default: 0)
      slope@double = 0;                 % Phase Slope (default: 0)
       
      FS@double = 4e7;                  % Sampling Frequency
      nFFT@uint16 = uint16(4096);       % Length of FFT (default: 4096)
      nSamples@uint16 = uint16(2030);   % Length of Time Domain Signal (default: 2030)
    end
    
    properties (Dependent = true)
       tdom@double; % Time Domain Impulse Response (length defined by nSamples)
       fdom@double; % Frequency Response (complex, length defined by nFFT)
    end
    
    methods
        % *****************************************************************
       %% CONSTRUCTOR
        function obj = msotEIR(varargin)
            % Create MSOT EIR Object
            %
            % msotEIR() - create with default properties
            % msotEIR(f_center,bw)
            % msotEIR(f_center,bw,phase_shift)
            % msotEIR(f_center,bw,phase_shift,slope)
            obj@handle;         % superclass constructor
            if nargin == 0, return; end
            
            if nargin >= 2,
                obj.f_center = varargin{1};
                obj.bw = varargin{2};
            end
            if nargin >= 3,
                obj.phase_shift = varargin{3};
            end
            if nargin >= 4,
                obj.slope = varargin{4};
            end
        end
        
        
        
        
        
        
        
        % *****************************************************************
       %% ACCESSORS
        function imp_resp = get.tdom(obj)
            imp_resp = calculateTD(obj);
        end
        
        function fdom = get.fdom(obj)
            fdom = calculateFD(obj);
        end
        
        function set.nFFT(obj,n)
            if mod(log2(double(n)),1) ~= 0,
                error('MSOT:msotEIR:InvalidNFFT','Length of FFT must be a power of 2');
            end
            obj.nFFT = n;
        end
        
        function str = toString(obj)
            str = sprintf('%.1fMHz @ %i%%, %.2fpi',...
                obj.f_center*1e-6,...
                obj.bw * 100,...
                obj.phase_shift);
        end
        
        
        
        
        % *****************************************************************
       %% CALCULATIONS
        
        function fdom = calculateFD(obj,varargin)
            % calculateFD  Calculate Frequency Domain Impulse Response
            % The Electric Impulse Response of the Transducer is calculated
            % from its Specifications. Per Default, the objects properties
            % are used for calculation - they can be overridden
            % using the optional parameters
            %
            % Parameters:
            % 1: (optional): Center Frequency in Hz 
            % 2: (optional): Relative T/R Bandwidth (0.6 = 60%)
            % 3: (optional): Phase Shift (default: 0)
            f = obj.f_center;
            b = obj.bw;
            p = obj.phase_shift;
            if nargin >= 3, f = varargin{1}; b = varargin{2}; end
            if nargin >= 4, p = varargin{3}; end
            
            % FWHM in frequency domain 
            fwhm = b*sqrt(2)*f;

            % gaussian width
            w=fwhm/(2*sqrt(2*log(2)));

            % frequency vector
            freq=obj.FS*linspace(0,1,obj.nFFT);

            % gaussian ampitude spectrum 
            mag_spec=exp(-((freq-f).^2/(2*w^2)));
            mag_spec=mag_spec+flip(mag_spec);

            % linear phase
            phase_spec=(linspace(0,1,obj.nFFT)*2*pi)*obj.slope+p;
%             phase_spec(uint16(obj.nFFT/2)+1:obj.nFFT) = -flip(phase_spec(uint16(obj.nFFT/2)+1:obj.nFFT));
            
%             figure;plot(freq,mag_spec);
%             figure;plot(freq,phase_spec);
            fdom = mag_spec.*exp(1i*phase_spec);
        end
        
        function imp_resp = calculateTD(obj,varargin)
            % calculateTD  Calculate Time Domain Impulse Response
            % The Electric Impulse Response of the Transducer is calculated
            % from its Specifications. 
            %
            % See also calculateFD
            fdom = calculateFD(obj,varargin{:});

            % transform to time domain
            iimp=real(ifft(fdom));

            % shift to center of time vector
            iimp=circshift(iimp',-int32(obj.nFFT / 2)+1);
            off = int32( ( numel(iimp) - obj.nSamples ) / 2);
            iimp = iimp(off+1:numel(iimp)-off);
%             figure;plot(iimp);

            % normalize
            norm_iimp=sqrt(sum(iimp.^2));
            imp_resp=(iimp/norm_iimp)';            
        end
        
        
        
        % *** Deconvolve the input Signal
        function sig = deconvolve(obj,sig,varargin)
            % deconvolve  Deconvolve Impulse Response From Signal
            % Deconvolution of Impulse Response from the Signal in
            % Frequency Domain. The length of the FFT is defined by nFFT.
            % 
            % Parameters:
            % (1): Signal (samples x channels)
            if isempty(sig), error('MSOT:signalProcessor:NoInputSignal','No input signal given'); end

            % waitbar
            wbar = [];
            if numel(varargin) >= 1,
                if isa(varargin{1},'logical') && varargin{1}, wbar = waitbar(0,'Filtering...'); 
                elseif isa(varargin{1},'matlab.ui.Figure'), wbar = varargin{1}; end
            end

            % dimension of input signal
            svec = size(sig);
            nsamples = svec(1);
            nchannels = svec(2);
            nsig = prod(svec(3:end));
            
            % deconvolve
            H = obj.fdom;
            sqrte = sqrt(eps);
            denom = abs(H).^2;
            ind = denom < sqrte;
            denom(ind)=sqrte;

            H = conj(H)./denom;
            HH = repmat(H', 1, nchannels);

            % transform to frequency domain if necessary
            if isreal(sig),
                if size(sig,1) > obj.nFFT,
                    warning('Length of Signal is longer than specified length of FFT');
                end
                if ~isempty(wbar), waitbar(0,wbar,'FFT...'); end;
                sig = sig - repmat(mean(sig,1),nsamples,1);
                sig = fft(sig,obj.nFFT,1);
                tdom = true;
            else
                tdom = false;
            end
            
            % apply filter
            sig = reshape(sig,[size(HH) nsig]);
            for jj = 1:nsig,
                if ~isempty(wbar), waitbar(jj/nsig,wbar,'Deconvolving...'); end;
                sig(:,:,jj) = sig(:,:,jj) .* HH;
            end;
            
            % transform back to time domain
            if tdom,
                if ~isempty(wbar), waitbar(0,wbar,'iFFT...'); end;
                sig = real(ifft(sig,[],1));
                sig = sig(1:nsamples,:);
            end
            
            sig = reshape(sig,svec);
            if numel(varargin) >= 1 && isa(varargin{1},'logical') && varargin{1}, close(wbar); end
        end

        
        % *****************************************************************
       %% CONFIGURATION
       function configure(obj)
           ret = inputdlg({'Center Frequency (MHz)','Bandwidth (% T/R)','Phase Offset'},'Impulse Response',1,...
                { num2str(obj.f_center*1e-6,'%.1f'),...
                  num2str(round(obj.bw*1e2),'%i'),...
                  num2str(obj.phase_shift,'%.2f')   });
           if isempty(ret), return; end;
           fc = str2double(ret{1}); if isnan(fc) || ~(fc > 1), errordlg('Cannot convert Frequency %s to a number (%.1f), aborting',ret{1},fc); return; end;
           baw = str2double(ret{2}); if isnan(baw) ||~(baw > 1), errordlg('Cannot convert bandwidth %s to a number (%.1f), aborting',ret{2},baw); return; end;
           ph = str2double(ret{3}); if isnan(ph), errordlg('Cannot convert phase %s to a number (%.2f), aborting',ret{3},ph); return; end;
           obj.f_center = fc*1e6;
           obj.bw = baw*1e-2;
           obj.phase_shift = ph;
       end
        
        
        
        % *****************************************************************
       %% FILE ACCESS
        
        function filename = write(obj,filename)
            % write   Write Time Domain IR to File
            % Writes Impulse Response to file compatible with ViewMSOT and
            % ViewMSOTc. If no parameters are given, a generic file name
            % based on the properties is used and returned
            %
            % Parameters:
            % (1): Filename [optional]

            if nargin == 1 || isempty(filename),
                filename = ['IR_' strrep(sprintf('%.1f',obj.f_center*1e-6),'.','c') 'MHz_'...
                            num2str(round(obj.bw*100)) 'perc' '.irf'];
            end

            FID = fopen(filename,'w');
            fwrite(FID,obj.tdom,'double');
            fclose(FID);
            fprintf('Impulse response written to %s\n',filename);
        end
        
    end
end