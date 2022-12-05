function [nlevel nlevelAvg m_noiseFilt m_noiseAvgFilt] = get_noiselevel(sigMat,navg,plot_dist,chplot,varargin)
% [nlevel nlevelAvg m_noiseFilt m_noiseAvgFilt] = get_noiselevel(sigMat,navg,plot_dist)
% 
% Inputs:
% * sigMat: signal matrix, dimensions: (2030,num channels,num samples)
% * navg: number of averages to evaluate (500 samples, 50 averages means 10
% sets to be evaluated
% * plot_dist: 1 to plot noise distribution graph, 0 otherwise
%
% Outputs:
% * nlevel: Mean noise level (single shot) in µV
% * nlevelAvg: Mean noise level (averaged) in µV
% * m_noiseFilt: Noise distribution (single shot) over channels in µV
% * m_noiseFiltAvg: Noise distribution (averaged) over channels in µV

%% parameters
setname = '';
if numel(varargin) >= 1
    setname = [varargin{1} ': '];
end
if numel(varargin) >= 2
    filter_f = varargin{2};
else
    filter_f = [200 10000] * 1e3;     % lower and higher cut-off frequencies in kHz
end



fprintf('\n******* Evaluating Signal noise level *******\n');
% sigrange = 700:size(sigMat,1);
sigrange = 700:1500;

%% *** PARAMETERS
if (min(sigMat(:)) < 1000 || max(sigMat(:)) > 4096)
    fprintf('Signal scaling is wrong, please use native DAQ range. Actual scaling: (%.0f - %0.f)\n',min(sigMat(:)),max(sigMat(:)));
%     return;
end
if (plot_dist)
    if (~isempty(setname))
        f = figure('Name',strrep(setname,':',''));
    else
        f = figure;
    end
    set(f,'Position',[200 100 800 1000]);
else f = [];
end

%% *** AVERAGE
% fprintf('* Averaging in sets of %i\n',navg);
% average different samples (non coherent noise)
for i = 1:floor(size(sigMat,3)/navg)
    sigMatAvg(:,:,i) = mean(sigMat(:,:,(i-1)*navg+1:i*navg),3);
end
clear i;

% *** average segments 
% TODO: Determine number of segments...
nseg = 4; 
if (size(sigMat,2) == 256)
    nseg = 8;
elseif (size(sigMat,2) == 384)
       nseg = 12;
elseif (size(sigMat,2) == 512)
       nseg = 16;
end
perseg = size(sigMat,2)/nseg;
sigMatSAvg = zeros(size(sigMat,1),nseg,size(sigMat,3));
for j = 1:nseg
    segstart(j) = (j-1)*perseg+1;
    segend(j) = j*perseg;
    sigMatSAvg(:,j,:) = mean(sigMat(:,segstart(j):segend(j),:),2);
end

%% *** FILTER 

% Filter setup
% fprintf('* Filtering signals (%.0f - %.0fHz)\n',filter_f);

fs = 40e6;                      % 40 Msamples DAQ sampling rate
len = size(sigMat, 1) ;
ts = 1/fs:1/fs:len/fs ;
dt = ts(2)-ts(1);

f_HPF = filter_f(1) ;
f_LPF = filter_f(2) ;
[b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .9 ) ;
% new filter
% [b_HPF,a_HPF] = cheby1( 1, .1, 2 * f_HPF/fs , 'high' ) ;
% old filter
[b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;
clear f_HPF f_LPF len ts dt;

% Filter all averaged signals
for i = 1:size(sigMatAvg,3)
    sigMatAvgFilt(:,:,i) = sigMatAvg(:,:,i) - ones(size(sigMatAvg,1),1)*mean(sigMatAvg(sigrange,:,i));
    sigMatAvgFilt(:,:,i) = FiltFiltM( b_LPF, a_LPF, sigMatAvgFilt(:,:,i), 1 ) ;
    sigMatAvgFilt(:,:,i) = FiltFiltM( b_HPF, a_HPF, sigMatAvgFilt(:,:,i), 1 ) ;
end

% Filter all signals
for i = 1:size(sigMat,3)
    sigMatFilt(:,:,i) = sigMat(:,:,i) - ones(size(sigMat,1),1)*mean(sigMat(sigrange,:,i));
    sigMatFilt(:,:,i) = FiltFiltM( b_LPF, a_LPF, sigMatFilt(:,:,i), 1 ) ;
    sigMatFilt(:,:,i) = FiltFiltM( b_HPF, a_HPF, sigMatFilt(:,:,i), 1 ) ;
end
% Filter all signals averaged across segments
for i = 1:size(sigMat,3)
    sigMatSAvgFilt(:,:,i) = sigMatSAvg(:,:,i) - ones(size(sigMatSAvg,1),1)*mean(sigMatSAvg(sigrange,:,i));
    sigMatSAvgFilt(:,:,i) = FiltFiltM( b_LPF, a_LPF, sigMatSAvgFilt(:,:,i), 1 ) ;
    sigMatSAvgFilt(:,:,i) = FiltFiltM( b_HPF, a_HPF, sigMatSAvgFilt(:,:,i), 1 ) ;
end
clear i b_LPF a_LPF b_HPF a_HPF filter_f

%% *** non-coherent NOISE LEVEL (averaged shots)

for i = 1:size(sigMatAvgFilt,3)
    noiseAvgFilt(i,:) = std(sigMatAvgFilt(sigrange,:,i)*4);
end
m_noiseAvgFilt = mean(noiseAvgFilt,1);
for j = 1:nseg
    nlevelAvg(j) = max(m_noiseAvgFilt(segstart(j):segend(j)));
end
fprintf(['Average Non-Coherent Noise level (' num2str(navg) ' avg, filtered): \t\t%.2fµV\n'],nlevelAvg);



%% *** coherent NOISE LEVEL on segment (averaged shots)

for i = 1:size(sigMatSAvgFilt,3)
    noiseSAvgFilt(i,:) = std(sigMatSAvgFilt(sigrange,:,i)*4);
end
nlevelSAvg = mean(noiseSAvgFilt,1);
fprintf(['Average Coherent Noise level (' num2str(navg) ' avg, filtered): \t\t%.2fµV\n'],nlevelSAvg);

%% *** Max NOISE LEVEL (single shot)

for i = 1:size(sigMatFilt,3)
    noiseFilt(i,:) = std(sigMatFilt(sigrange,:,i)*4);
end
m_noiseFilt = mean(noiseFilt,1);
for j = 1:nseg
    nlevel(j) = mean(m_noiseFilt(segstart(j):segend(j)));
end
fprintf('Max Noise level (single shot, filtered): \t%.2fµV\n',nlevel);



%% graph...
if (~isempty(f))

    % non-coherent noise
    subplot(2,1,1);
    semilogy(m_noiseAvgFilt,'Color','black','LineWidth',1,'DisplayName',['A: non-coherent noise across ' num2str(navg) ' samples']);
    hold on;

    % coherent segment noise
    xplot = [];
    yplot = [];
    for j = 1:nseg
        xplot = [xplot segstart(j):segend(j)];
        yplot = [yplot; ones(perseg,1).*nlevelSAvg(j)];
    end
    plot(xplot,yplot,'Color',[0 0.5 0],'DisplayName','B: coherent segment noise');
    
    
    % single shot noise
    semilogy(m_noiseFilt,'Color','blue','LineWidth',1,'DisplayName',['C: single shot (filtered)']);

    
    legend('show');
    grid on;
    set(gca,'XTick',segend)
    % legend('show','location','EastOutside');
%     line([1 size(sigMat,2)],[nlevel nlevel],'Color','blue','LineStyle',':');
%     text(size(sigMat,2)-2,nlevel,['mean: ' num2str(nlevel,'%.2f') 'µV'],'Color','blue','BackgroundColor','white','HorizontalAlign','right');
%     line([1 size(sigMat,2)],[nlevelAvg nlevelAvg],'Color','black','LineStyle',':');
%     text(size(sigMat,2)-2,nlevelAvg,['mean: ' num2str(nlevelAvg,'%.2f') 'µV'],'Color','black','BackgroundColor','white','HorizontalAlign','right');
    set(gca,'ycolor',[0.4 0.4 0.4]);
    ylabel('Standard deviation [µV]');
    xlabel('channel');
    ylim([0.7 80]);
    xlim([1 size(sigMat,2)]);
    
    line(xlim,[1 1]*2,'Color','black','LineStyle','--','LineWidth',1);
    line(xlim,[1 1]*6,'Color',[0 0.5 0],'LineStyle','--','LineWidth',1);
    line(xlim,[1 1]*12,'Color','blue','LineStyle','--','LineWidth',1);
    set(gca,'XColor',[0.3 0.3 0.3]);
    set(gca,'YColor',[0.7 0.7 0.7]);
    text(2,2,'2µV','Color','black','VerticalAlign','bottom');
    text(2,6,'6µV','Color',[0 0.5 0],'VerticalAlign','bottom');
    text(2,12,'12µV','Color','blue','VerticalAlign','bottom');
    title([setname 'Noise distribution (' num2str(size(sigMat,3)) ' samples)'],'Interpreter','none');
    
    
    
    
    % Noise distribution
    subplot(2,1,2);
    imagesc(noiseFilt);
    ch=colorbar;
    xlabel('channels');
    ylabel('samples');
    title([setname 'standard deviation of noise sample']);
    axis image;
    set(get(ch,'YLabel'),'String','Noise STD [µV]');
    colormap(jet);

    caxis([8 12]);

end

%% plot individual channels
if (~isempty(chplot))
   g = figure;
    t = (1/fs:1/fs:size(sigMat,1)/fs)*1e6;
   % plot every channel
   for i = 1:numel(chplot)
       sig = sigMatAvgFilt(:,chplot(i),1);
       
       subplot(numel(chplot),2,(i-1)*2+1);
       plot(t,sig);
       title(['Signal channel ' num2str(chplot(i)) ' (' num2str(navg) ' averages)']);
       xlim([t(1) t(end)]);
       %ylim([mean(sig(:))-100 mean(sig(:))+100]);
       xlabel('t [µs]');
       ylabel('signal (a.u.)');

       subplot(numel(chplot),2,(i-1)*2+2);
       [Y f NFFT] = calc_fft(sig -(mean(sig)));
       semilogy(f*1e-6,2*abs(Y(1:NFFT/2+1))); 
       title('Single sided amplitude spectrum');
       xlim([0 10]);
       xlabel('frequency [MHz]');
       ylim([0.001 5]);
       ylabel('|Y(f)|')

   end 
end
