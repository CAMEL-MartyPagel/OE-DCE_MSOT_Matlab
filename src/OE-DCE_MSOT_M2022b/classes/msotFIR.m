classdef msotFIR < msotFilter
    % msotFIR   View MSOT Frequency Domain FIR Implementation
    % Class that mimics OpenCL FIR Implementation (created alonside the
    % msotFIR method supplied by Mike)
    properties 
        differentiator@logical = true;      % Add complex slope for differentiator
        fHP@double = double(50e3);          % Low Edge of Passband
        fLP@double = double(6e6);           % High Edge of Passband
        filter_rise@double = double(0.2);   % Filter Rise in %
        filter_fall@double = double(0.2);   % Filter Fall in %
        
        order@double = double(1022);        % FIR order
        npt@double = double(2048);          % number of Points for Interpolation of Spectrum
        conjugation@logical = true;         % Conjugation of the numerator

        fVec@double = [];
        aVec@double = [];
    end
    
    methods
        function fv = get.fVec(obj)
            % calculate frequency vector
            if isempty(obj.fVec),
                fv = [0, (1-obj.filter_rise)*obj.fHP, obj.fHP,     obj.fLP,     (1+obj.filter_fall)*obj.fLP, obj.FS]/obj.FS;
%                 error('MSOT:msotFIR:NoCharacteristics','No Frequency Vector or Characteristics set in msotFIR');                
            else
                fv = obj.fVec;
            end
        end
        
        function av = get.aVec(obj)
            % calculate frequency vector
            if isempty(obj.aVec),
                ampmin = 1; 
                ampmax = 1;
                if obj.differentiator,
                    ampmin = (1/(obj.FS)*obj.fHP)*1i;
                    ampmax = (1/(obj.FS)*obj.fLP)*1i;
                end                    
                av = double([0, 0, ampmin, ampmax, 0, 0]);
%                 error('MSOT:msotFIR:NoCharacteristics','No Frequency Vector or Characteristics set in msotFIR');                
            else
                av = obj.aVec;
            end
        end
        
        
        
        function str = toString(obj)
            le = sprintf('%iHz',obj.fHP);
            if obj.fHP > 1e3, le = sprintf('%ikHz',round(obj.fHP/1e3)); end
            if obj.fHP > 1e6, le = sprintf('%.1fMHz',obj.fHP/1e6); end
            
            he = sprintf('%iHz',obj.fLP);
            if obj.fLP > 1e3, he = sprintf('%ikHz',round(obj.fLP/1e3)); end
            if obj.fLP > 1e6, he = sprintf('%.1fMHz',obj.fLP/1e6); end

            str = sprintf('%s: %s-%s','FIR',le,he);
        end
        
        
        % *****************************************************************
       %% CONFIGURATION
       function configure(obj)
           ret = inputdlg({'Lower Cutoff (Hz)','Upper Cutoff (Hz)'},'FIR Filter',1,...
                { num2str(round(obj.fHP),'%i'),...
                  num2str(round(obj.fLP),'%i') });
           if isempty(ret), return; end;
           hpf = str2double(ret{1}); if isnan(hpf) , errordlg(sprintf('Cannot convert Frequency %s to a number (%.1f), aborting',ret{1},hpf)); return; end;
           lpf = str2double(ret{2}); if isnan(lpf) || ~(lpf > 0), errordlg(sprintf('Cannot convert bandwidth %s to a number (%.1f), aborting',ret{2},lpf)); return; end;
           obj.fHP = hpf;
           obj.fLP = lpf;
       end
        
        
        % *****************************************************************
       %% CALCULATION
        function fd = fdom(obj)
            % fdom   Calculate Frequency Domain Response according to Settings
            aa = obj.aVec;
            ff = obj.fVec;
            
            % validate Parameters
            nn = double(obj.order + 1);
            if ((abs(aa(end)) ~= 0) && (ff(end) == 1.0) && (mod(obj.order,2) == 1))
                nn = nn + 1;
            end

            df = diff(ff);
            chk = df<0;
            if sum(chk)>0,
                error('MSOT:msotFIR:InvalidFVec','interpolateSpectra: Invalid Frequency Vector');
            end

            
            lap = 16;
            sizeH = 2*obj.npt-2;
            H = zeros(sizeH,1);

            % interpolate Spectrm
            nb = 2;
            H(1) = aa(1);
            for i=1:numel(ff)-1,
                if df(i) == 0,
                    nb = nb - floor((lap-1)/2);
                    ne = nb + lap;
                else
                    ne = floor(ff(i+1)*double(obj.npt));
                end;
                if nb < 1 || ne > obj.npt,
                    error('MSOT:msotFIR:SpecError','interpolateSpectra: Specification Error');
                end;
                if nb == ne
                    inc = 0;
                else
                    j=nb:ne;
                    inc = (j-nb)/(ne-nb);
                end
                H(nb:ne) = inc*aa(i+1) + (1 - inc)*aa(i);
                nb = ne + 1;
            end

            % Fourier time-shift.
            k = double(1:obj.npt);
            fac = -0.5i*pi * (nn - 1) / (double(obj.npt)-1)*(k'-1);
            H(k) = H(k) .* exp(fac);

            k = 2:(obj.npt-1);
            if (obj.conjugation == true)
                rH = conj(H(k));
            else
                rH = H(k);
            end
            H(sizeH-k+2) = rH;
            % end

            % time domain response and hamming windowing thereof
            ht = real(ifft(H));
            win = hamming(nn);
            b = ht(int32(1:nn)) .* win;

            half = floor(nn/2);
            B = [b(half+1:end);zeros(4096-nn,1);b(1:half)];
            fd = fft(B);

        end
    end 
end