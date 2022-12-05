function rmsimg = get_stddist(sigMat,varargin)
% rmsimg = get_stddist(sigMat)
% visualize standard deviation across one noise sample.
% 
% Inputs:
% * sigMat: signal matrix, dimensions: (2030,num channels,num samples)
%
% Outputs:
% * rmsimg

%% parameters
setname = '';
if numel(varargin) >= 1
    setname = [varargin{1} ': '];
end

fprintf('\n******* Evaluating Max Signal Noise Sample *******\n');
% sigrange = 700:size(sigMat,1);
sigrange = 1:2030;

%% *** PARAMETERS
if (min(sigMat(:)) < 1000 || max(sigMat(:)) > 4096)
    fprintf('Signal scaling is wrong, please use native DAQ range. Actual scaling: (%.0f - %0.f)\n',min(sigMat(:)),max(sigMat(:)));
%     return;
end
f = figure;



%% *** FILTER 

% Filter setup
filter_f = [50 10000] * 1e3;     % lower and higher cut-off frequencies in kHz
% fprintf('* Filtering signals (%.0f - %.0fHz)\n',filter_f);

fs = 40e6;                      % 40 Msamples DAQ sampling rate
len = size(sigMat, 1) ;
ts = 1/fs:1/fs:len/fs ;
dt = ts(2)-ts(1);

f_HPF = filter_f(1) ;
f_LPF = filter_f(2) ;
[b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .9 ) ;
[b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;
clear f_HPF f_LPF len ts dt;


% Filter all signals
for i = 1:size(sigMat,3)
    sigMatFilt(:,:,i) = sigMat(:,:,i) - ones(size(sigMat,1),1)*mean(sigMat(sigrange,:,i));
    sigMatFilt(:,:,i) = FiltFiltM( b_LPF, a_LPF, sigMatFilt(:,:,i), 1 ) ;
    sigMatFilt(:,:,i) = FiltFiltM( b_HPF, a_HPF, sigMatFilt(:,:,i), 1 ) ;
end
clear i b_LPF a_LPF b_HPF a_HPF filter_f

%% *** Max NOISE LEVEL (single shot)

for i = 1:size(sigMatFilt,3)
    rmsimg(i,:) = std(sigMatFilt(sigrange,:,i)*4);
end

imagesc(rmsimg);
ch=colorbar;
xlabel('channels');
ylabel('samples');
title([setname 'standard deviation of noise sample']);
axis image;
set(get(ch,'YLabel'),'String','Noise STD [µV]');

fprintf('Minimum STD: \t%.1fµV\n',min(rmsimg(:)));
fprintf('Mean STD: \t\t%.1fµV\n',mean(rmsimg(:)));
fprintf('Maximum STD: \t%.1fµV\n',max(rmsimg(:)));
% caxis([min(rmsimg(:)) 12]);
caxis([8 12]);
