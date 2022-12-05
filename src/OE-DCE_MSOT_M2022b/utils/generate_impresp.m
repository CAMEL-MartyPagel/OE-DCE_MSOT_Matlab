function imp_resp = generate_impresp(f_center,bw,varargin)
% GENERATE_IMPRESP Generate Impulse Response
%    imp_resp = generate_impresp(f_center,bw,[phase],[slope],[filename])
%
% Parameters:
% - f_center: Center frequency in Hz, e.g. 5e6;
% - bw: Fractional bandwidth (send/receive) - e.g. 0.6 for 60%
% - (optional) phase: phase offset - e.g. 0.5 for 0.5*pi
% - (optional) slope: linear phase slope 
% - (optional) filename: Filename to write XVUE file, [] for auto-naming
% 
% usage: 
% imp_resp = generate_impresp(f_center,bw)
%   Generate impulse response with certain center frequency and bandwidth
% 
% imp_resp = generate_impresp(f_center,bw,phase,slope)
%   Specify phase offset in fractions of pi, and slope
% 
% imp_resp = generate_impresp(f_center,bw,phase,slope,filename)
%   Additionally write impulse response to file in XVUE format. Supply
%   emtpy array for filename to auto-create a name

% additional input parameters
if (numel(varargin) > 0 && ~isempty(varargin{1}))
    phase_shift=varargin{1}*pi;
else
    phase_shift=0;
end

if (numel(varargin) > 1 && ~isempty(varargin{2}))
    slope=varargin{2};
else
    slope=0;
end

% sampling rate
fs = 40e6;

% FWHM in frequency domain 
fwhm = bw*sqrt(2)*f_center;

% gaussian width
w=fwhm/(2*sqrt(2*log(2)));

% frequency vector
freq=fs*linspace(0,1,2030);

% gaussian ampitude spectrum 
mag_spec=exp(-((freq-f_center).^2/(2*w^2)));
mag_spec=mag_spec+flip(mag_spec);
% mag_spec(1:1015) = 0;

% linear phase
phase_spec=(linspace(0,1,2030)*2*pi)*slope+phase_shift;
phase_spec(1016:2030) = -flip(phase_spec(1016:2030));


% figure;plot(mag_spec);hold on;plot(phase_spec);
% transform to time domain
iimp=real(ifft(mag_spec.*exp(1i*phase_spec)));
% iimp_i=imag(ifft(mag_spec.*exp(1i*phase_spec)));
% iimp_a=abs(ifft(mag_spec.*exp(1i*phase_spec)));

% shift to center of time vector
iimp=circshift(iimp',-1015);
% iimp_i=circshift(iimp_i',-1015);
% iimp_a=circshift(iimp_a',-1015);
% figure;plot(iimp);hold on;plot(iimp_i);plot(iimp_a);
% title(num2str(phase_shift));xlim([1000 1030]);

% normalize
norm_iimp=sqrt(sum(iimp.^2));
imp_resp=(iimp/norm_iimp)';

% write to file
if (numel(varargin) > 2),
    if isnumeric(varargin{2}) && numel(varargin) > 2
        if isempty(varargin{3}),
            filename = ['IR_' strrep(sprintf('%.1f',f_center*1e-6),'.','c') 'MHz_'...
                        num2str(round(bw*100)) 'perc' '.irf'];
%             filename = ['IR_' strrep(sprintf('%.1f',f_center*1e-6),'.','c') 'MHz_'...
%                         num2str(round(bw*100)) 'perc_0d' num2str(round(phase_shift/pi*100),'%i') 'pi_slope' num2str(slope,'%i') '.irf'];
        else
            filename = varargin{3};
        end
    else
        if isempty(varargin{2}),
            filename = ['IR_' strrep(sprintf('%.1f',f_center*1e-6),'.','c') 'MHz_'...
                        num2str(round(bw*100)) 'perc' '.irf'];
%             filename = ['IR_' strrep(sprintf('%.1f',f_center*1e-6),'.','c') 'MHz_'...
%                         num2str(round(bw*100)) 'perc_0d' num2str(round(phase_shift/pi*100),'%i') 'pi_slope' num2str(slope,'%i') '.irf'];
        elseif ischar(varargin{2}),
            filename = varargin{2};
        else
            return;
        end
    end

    FID = fopen(filename,'w');
    fwrite(FID,imp_resp(1:2030)','double');
    fclose(FID);
    fprintf('Impulse response written to %s\n',filename);
end