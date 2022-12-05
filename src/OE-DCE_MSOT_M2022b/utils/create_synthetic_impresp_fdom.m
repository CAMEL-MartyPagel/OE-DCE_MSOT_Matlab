clear all;
%% calculate artificial impulse response III

% Specs according to Manufacturer
f_center=5.2e6;
bandwidth = 0.70;                       % send/receive
phase_shift=0.5*pi;

% Derived parameters
fwhm = bandwidth*sqrt(2)*f_center;       % FWHM of IR (receive only)
flip = 0;
fs = 40e6;
freq=fs*linspace(0,1,2030);
w=fwhm/(2*sqrt(2*log(2)));               % STD

% spectrum
mag_spec=exp(-((freq-f_center).^2/(2*w^2)));
mag_spec=mag_spec+flipud(mag_spec);
phase_spec=(zeros(1,2030))+phase_shift;

% impulse response
iimp=real(ifft(mag_spec.*exp(1i*phase_spec)));
iimp=circshift(iimp',-1015);

norm_iimp=sqrt(sum(iimp.^2));
imp_resp=(iimp/norm_iimp)';

fsig=fft(imp_resp);

irname = ['IR_synthetic_phi0_fdom_' strrep(num2str(f_center*1e-6,'%.1f'),'.','c') 'MHz_' num2str(bandwidth*100,'%i') 'perc_0c' num2str(round(phase_shift/pi*100),'%i') 'pi'];
FID = fopen([irname '.irf'],'w');
fwrite(FID,imp_resp(1:2030)'./sum(abs(imp_resp(1:2030))),'double');
fclose(FID);
