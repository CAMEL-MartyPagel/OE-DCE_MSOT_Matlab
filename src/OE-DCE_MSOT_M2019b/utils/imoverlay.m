function cimg = imoverlay(bgimg,fgimg,cpar);

par.fgmap = 'green';
par.fgathres = [0 0.75];
par.fgfilt = [0 0];
par.fgthres = [];
par.fgimg2 = [];
par.fgthres2 = [];
par.fgathres2 = [0 0.75];
par.fgmap2 = 'red';
par.bgthres = [];
par.bgfilt = [0 0];
par.yticks = [];
par.xticks = [];
par.sensimg = [];
par.tickalpha = 0.5;

fx = fieldnames(cpar);
for j = 1:numel(fx)
    par = setfield(par,fx{j},getfield(cpar,fx{j}));
end
clear cpar j fx;

fgathres = par.fgathres;
fgfilt = par.fgfilt;
fgthres = par.fgthres;
bgfilt = par.bgfilt;
bgthres = par.bgthres;
fgathres = par.fgathres;

n = size(fgimg,1);

% colormap
fgmap = getcmap(par.fgmap);

% generate background image
% bgimg = filterImage(bgimg,10.^bgfilt);

if (isempty(bgthres)), bgthres = autothres(bgimg); end;
bgimg = bgimg - bgthres(1); bgimg (bgimg < 0) = 0;
bgimg = bgimg ./ (bgthres(2)-bgthres(1)); bgimg(bgimg > 1) = 1;
% bgimg = 10*log10(bgimg);
% bgimg = bgimg - min(bgimg(:));
% bgimg = bgimg ./ max(bgimg(:));
bgimg_c = ind2rgb(round(bgimg*64),gray);




% fgimg = filterImage(fgimg,10.^fgfilt);
fgaimg = fgimg;     % save for alpha thresholding

% threshold and alpha threshold
if (isempty(fgthres)), fgthres = autothres(fgimg); end;
fgimg = fgimg - fgthres(1); fgimg (fgimg < 0) = 0;
fgimg = fgimg ./ (fgthres(2)-fgthres(1)); fgimg(fgimg > 1) = 1;

fgimg_c = ind2rgb(round(fgimg*64),fgmap);

% melanin alpha
fgathres = [min(fgthres)+diff(fgthres)*fgathres(1) min(fgthres)+diff(fgthres)*fgathres(2)];
fgaimg = fgaimg - fgathres(1); fgaimg (fgaimg < 0) = 0;
fgaimg = fgaimg ./ (fgathres(2)-fgthres(1)); fgaimg(fgaimg > 1) = 1;    

if ~isempty(par.fgimg2)
    fgmap2 = getcmap(par.fgmap2);
    fgimg2 = par.fgimg2;
    fgaimg2 = fgimg2;     % save for alpha thresholding
    fgthres2 = par.fgthres2;
    fgathres2 = par.fgathres2;

    % threshold and alpha threshold
    if (isempty(fgthres2)), fgthres2 = autothres(fgimg2); end;
    fgimg2 = fgimg2 - fgthres2(1); fgimg2 (fgimg2 < 0) = 0;
    fgimg2 = fgimg2 ./ (fgthres2(2)-fgthres2(1)); fgimg2(fgimg2 > 1) = 1;

    fgimg2_c = ind2rgb(round(fgimg2*64),fgmap2);

    % melanin alpha
    fgathres2 = [min(fgthres2)+diff(fgthres2)*fgathres2(1) min(fgthres2)+diff(fgthres2)*fgathres2(2)];
    fgaimg2 = fgaimg2 - fgathres2(1); fgaimg2 (fgaimg2 < 0) = 0;
    fgaimg2 = fgaimg2 ./ fgathres2(2); fgaimg2(fgaimg2 > 1) = 1;    
    
    % mixing
%     fg1 = reshape(fgimg_c,size(fgimg,1)*size(fgimg,2)*size(fgimg_c,3),1);
%     fg2 = reshape(fgimg2_c,size(fgimg,1)*size(fgimg,2)*size(fgimg_c,3),1);
%     fg = max([fg1 fg2],[],2);
%     fg = reshape(fg,size(fgimg_c));
%     
%     fga1 = reshape(fgaimg,size(fgaimg,1)*size(fgaimg,2),1);
%     fga2 = reshape(fgaimg2,size(fgaimg,1)*size(fgaimg,2),1);
%     fga = sum([fga1 fga2],2);
%     fga(fga > 1) = 1;
%     fga = reshape(fga,size(fgaimg));
%     
%     fgimg_c = fg;
%     fgaimg = fga;
    
%     fgaimg2 = fgaimg2 * 0.5;
end


f = figure('Visible','on','Position',[50 50 size(fgimg,2)+20 size(fgimg,1)+20]);
ax=axes('Units','pixels','Position',[10 10 size(fgimg,2) size(fgimg,1)]);

image(bgimg_c); 
hold on;
axis off;


% first overlay
ih = image(fgimg_c);
set(ih,'AlphaData',fgaimg);

if ~isempty(par.fgimg2)
    ih = image(fgimg2_c);
    set(ih,'AlphaData',fgaimg2);
end

tickimg = zeros(size(fgimg));
if ~isempty(par.sensimg),
    tickimg = tickimg + bwperim(logical(par.sensimg));
end

for y = par.yticks,
    tickimg(y,[1:10 end-10:end]) = 1;
%     tickimg(y,[1:4:end 2:4:end]) = 1;
end

for x = par.xticks,
    tickimg([1:10 end-10:end],x,1) = 1;
%     tickimg([1:4:end 2:4:end],x,1) = 1;
end

if nnz(tickimg(:)) > 0
    ih1 = image(ind2rgb(tickimg*64,gray));
    set(ih1,'AlphaData',tickimg*par.tickalpha);
end

% Copy Frame
fr = getframe(ax);
cimg = fr.cdata(2:2+n-1,1:1+n-1,:);
% cimg=hardcopy(gca,'-Dzbuffer', '-r0');
% cimg = cimg(end-10-n:end-10-1,13:n+12,:);
close(f);



% figure;image(cimg);

return;
xt = linspace(0,size(fgimg,2),5);
xt(1) = 1;
set(gca,'XTick',xt);
set(gca,'XTickLabel',num2str(linspace(-size(fgimg,2)/2*0.2,size(fgimg,2)/2*0.2,5)'));
set(gca,'YTick',xt);
set(gca,'YTickLabel',num2str(linspace(-size(fgimg,2)/2*0.2,size(fgimg,2)/2*0.2,5)'));
xlabel('x [mm]');
ylabel('y [mm]');
set(gca,'XColor',[0.4 0.4 0.4]);
set(gca,'YColor',[0.4 0.4 0.4]);

line([20 45],[size(fgimg,1)-10 size(fgimg,1)-10],'Color','white','LineWidth',2);
text(33,size(fgimg,1)-20,'5mm','Color','white','HorizontalAlign','center','FontSize',8);