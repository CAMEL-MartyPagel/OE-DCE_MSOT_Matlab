function cimg = imoverlay3(bgimg,fgimg,cpar,varargin)

par.fgmap = 'green';
par.fgathres = [0 0.75];
par.fgfilt = [0 0];
par.fgthres = [];
par.fgimg2 = [];
par.fgthres2 = [];
par.fgathres2 = [0 0.75];
par.fgmap2 = 'green';
par.bgthres = [];
par.bgfilt = [0 0];
par.yticks = [];
par.xticks = [];
par.zticks = [];
par.pos = [0 0 0];      % [0 0 0] for MIP, otherwise positions
par.scale = 1;
par.tickalpha = 0.5;
par.sensvol = [];
par.senscolor = [100 140 200]./255;
par.aspect = [1 1 1];
par.scalebar = [];
par.roi = [];
par.roicolor = [120 230 235]./255;

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


% size(par.roi)
% interpolation
if prod(par.aspect) ~= 1,
    [x0 y0 z0] = meshgrid(linspace(0,1,size(fgimg,1)),...
        linspace(0,1,size(fgimg,2)),...
        linspace(0,1,size(fgimg,3)));
    [x_ y_ z_] = meshgrid(linspace(0,1,size(fgimg,1)*par.aspect(1)),...
        linspace(0,1,size(fgimg,2)*par.aspect(2)),...
        linspace(0,1,size(fgimg,3)*par.aspect(3)));
    bgimg = interp3(x0,y0,z0,bgimg,x_,y_,z_);
    fgimg = interp3(x0,y0,z0,fgimg,x_,y_,z_);
    if ~isempty(par.sensvol),
        par.sensvol = logical(interp3(x0,y0,z0,double(par.sensvol),x_,y_,z_));
    end
    if ~isempty(par.roi),
        for jj = 1:size(par.roi,4),
            roi2(:,:,:,jj) = logical(interp3(x0,y0,z0,double(par.roi(:,:,:,jj)),x_,y_,z_));
        end
        par.roi = roi2; 
        clear roi2;
    end;
    
    if ~isempty(par.fgimg2),
        par.fgimg2 = interp3(x0,y0,z0,par.fgimg2,x_,y_,z_);
    end
    clear x0 y0 z0 x_ y_ z_;
end

n = size(fgimg,1);
nz = size(fgimg,3);

spac = 10;

% colormaps
fgmap = getcmap(par.fgmap);


% Thresholding
if (isempty(bgthres)), bgthres = autothres(bgimg); end;
bgimg = bgimg - bgthres(1); bgimg (bgimg < 0.0) = 0.0;
bgimg = bgimg ./ (bgthres(2) - bgthres(1)); bgimg(bgimg > 1) = 1;
% bgimg = 10*log10(bgimg);
% bgimg = bgimg - min(bgimg(:));
% bgimg = bgimg ./ max(bgimg(:));


fgaimg = fgimg;     % save for alpha thresholding

% threshold and alpha threshold
if (isempty(fgthres)), fgthres = autothres(fgimg); end;
fgimg = fgimg - fgthres(1); fgimg (fgimg < 0) = 0;
fgimg = fgimg ./ (fgthres(2) - fgthres(1)); fgimg(fgimg > 1) = 1;

fgathres = [min(fgthres)+diff(fgthres)*fgathres(1) min(fgthres)+diff(fgthres)*fgathres(2)];
fgaimg = fgaimg - fgathres(1); fgaimg (fgaimg < 0) = 0;
fgaimg = fgaimg ./ (fgathres(2)- fgathres(1)); fgaimg(fgaimg > 1) = 1;    


if ~isempty(par.fgimg2),
    fgimg2 = par.fgimg2;
    fgaimg2 = fgimg2;
    fgmap2 = getcmap(par.fgmap2);
    
    % threshold and alpha threshold
    if (isempty(par.fgthres2)), par.fgthres2 = autothres(fgimg2); end;
    fgthres2 = par.fgthres2;
    fgimg2 = fgimg2 - par.fgthres2(1); fgimg2 (fgimg2 < 0) = 0;
    fgimg2 = fgimg2 ./ par.fgthres2(2); fgimg2(fgimg2 > 1) = 1;

    fgathres2 = [min(fgthres2)+diff(fgthres2)*par.fgathres2(1) min(fgthres2)+diff(fgthres2)*par.fgathres2(2)];
    fgaimg2 = fgaimg2 - fgathres2(1); fgaimg2 (fgaimg2 < 0) = 0;
    fgaimg2 = fgaimg2 ./ (fgathres2(2)-fgathres2(1)); fgaimg2(fgaimg2 > 1) = 1;    
    
    
end








% cross-section or projection
if prod(par.pos) > 0,
    bgxy = bgimg(:,:,par.pos(3));
    fgxy = fgimg(:,:,par.pos(3));
    fgaxy = fgaimg(:,:,par.pos(3));
    bgzx = squeeze(bgimg(par.pos(1),:,:));
    fgzx = squeeze(fgimg(par.pos(1),:,:));
    fgazx = squeeze(fgaimg(par.pos(1),:,:));
    bgyz = squeeze(bgimg(:,par.pos(2),:))';
    fgyz = squeeze(fgimg(:,par.pos(2),:))';
    fgayz = squeeze(fgaimg(:,par.pos(2),:))';
else
    bgxy = max(bgimg,[],3);
    fgxy = max(fgimg,[],3);
    fgaxy = max(fgaimg,[],3);
    bgyz = squeeze(max(bgimg,[],1))';
    fgyz = squeeze(max(fgimg,[],1))';
    fgayz = squeeze(max(fgaimg,[],1))';
    bgzx = squeeze(max(bgimg,[],2));
    fgzx = squeeze(max(fgimg,[],2));
    fgazx = squeeze(max(fgaimg,[],2));

    if ~isempty(par.fgimg2),
        fgxy2 = max(fgimg2,[],3);
        fgaxy2 = max(fgaimg2,[],3);
        fgyz2 = squeeze(max(fgimg2,[],1))';
        fgayz2 = squeeze(max(fgaimg2,[],1))';
        fgzx2 = squeeze(max(fgimg2,[],2));
        fgazx2 = squeeze(max(fgaimg2,[],2));
    end
    
end

bgxy_c = ind2rgb(round(bgxy*64),gray);
bgyz_c = ind2rgb(round(bgyz*64),gray);
bgzx_c = ind2rgb(round(bgzx*64),gray);
fgxy_c = ind2rgb(round(fgxy*64),fgmap);
fgyz_c = ind2rgb(round(fgyz*64),fgmap);
fgzx_c = ind2rgb(round(fgzx*64),fgmap);
if ~isempty(par.fgimg2),
    fgxy2_c = ind2rgb(round(fgxy2*64),fgmap2);
    fgyz2_c = ind2rgb(round(fgyz2*64),fgmap2);
    fgzx2_c = ind2rgb(round(fgzx2*64),fgmap2);
    
    % Color Mix the two layers
    fgxy_c = reshape(max([...
        reshape(fgxy_c,size(fgxy,1)*size(fgxy,2)*size(fgxy_c,3),1), ...
        reshape(fgxy2_c,size(fgxy,1)*size(fgxy,2)*size(fgxy_c,3),1)...
        ],[],2),size(fgxy_c));
    fgaxy = reshape(max([...
            reshape(fgaxy,size(fgaxy,1)*size(fgaxy,2),1),...
            reshape(fgaxy2,size(fgaxy,1)*size(fgaxy,2),1)...
        ],[],2),size(fgaxy));
    
    fgyz_c = reshape(max([...
        reshape(fgyz_c,size(fgyz,1)*size(fgyz,2)*size(fgyz_c,3),1), ...
        reshape(fgyz2_c,size(fgyz,1)*size(fgyz,2)*size(fgyz_c,3),1)...
        ],[],2),size(fgyz_c));
    fgayz = reshape(max([...
            reshape(fgayz,size(fgayz,1)*size(fgayz,2),1),...
            reshape(fgayz2,size(fgayz,1)*size(fgayz,2),1)...
        ],[],2),size(fgayz));
     
    fgzx_c = reshape(max([...
        reshape(fgzx_c,size(fgzx,1)*size(fgzx,2)*size(fgzx_c,3),1), ...
        reshape(fgzx2_c,size(fgzx,1)*size(fgzx,2)*size(fgzx_c,3),1)...
        ],[],2),size(fgzx_c));
    fgazx = reshape(max([...
            reshape(fgazx,size(fgazx,1)*size(fgazx,2),1),...
            reshape(fgazx2,size(fgazx,1)*size(fgazx,2),1)...
        ],[],2),size(fgazx));
     
    
    
end

% Add dotted lines
for y = par.yticks,
    fgxy_c(y,[1:10 end-10:end],:) = 255;
    fgaxy(y,[1:10 end-10:end],:) = par.tickalpha;
    fgyz_c([1:10 end-10:end],y,:) = 255;
    fgayz([1:10 end-10:end],y,:) = par.tickalpha;
end

for x = par.xticks,
    fgxy_c([1:10 end-10:end],x,:) = 255;
    fgaxy([1:10 end-10:end],x,:) = par.tickalpha;
    fgzx_c(x,[1:10 end-10:end],:) = 255;
    fgazx(x,[1:10 end-10:end],:) = par.tickalpha;
end

for z = par.zticks,
    fgzx_c([1:10 end-10:end],z,:) = 255;
    fgazx([1:10 end-10:end],z,:) = par.tickalpha;
    fgyz_c(z,[1:10 end-10:end],:) = 255;
    fgayz(z,[1:10 end-10:end],:) = par.tickalpha;
end

if ~isempty(par.scalebar),
    fgxy_c([end-13:end-10],[end-10-par.scalebar:end-10],:) = 255;
    fgaxy([end-13:end-10],[end-10-par.scalebar:end-10],:) = 1;
end
% line([n-10-scalebarwidth n-10],[n-10 n-10],'LineWidth',4,'Color','white');
        

if ~isempty(par.sensvol),
    sensxy = max(logical(par.sensvol),[],3);
    sensxy = bwperim(sensxy); 
    i = fgxy_c(:,:,1); i(sensxy) = par.senscolor(1);fgxy_c(:,:,1) = i;
    i = fgxy_c(:,:,2); i(sensxy) = par.senscolor(2);fgxy_c(:,:,2) = i;
    i = fgxy_c(:,:,3); i(sensxy) = par.senscolor(3);fgxy_c(:,:,3) = i;
    fgaxy(sensxy) = 1;
    
    sensyz = squeeze(max(logical(par.sensvol),[],1))';
    sensyz = bwperim(sensyz);
    i = fgyz_c(:,:,1); i(sensyz) = par.senscolor(1);fgyz_c(:,:,1) = i;
    i = fgyz_c(:,:,2); i(sensyz) = par.senscolor(2);fgyz_c(:,:,2) = i;
    i = fgyz_c(:,:,3); i(sensyz) = par.senscolor(3);fgyz_c(:,:,3) = i;
    fgayz(sensyz) = 1;

    senszx = squeeze(max(logical(par.sensvol),[],2));
    senszx = bwperim(senszx);[x,y]=ind2sub(size(senszx), find(senszx));
    i = fgzx_c(:,:,1); i(senszx) = par.senscolor(1);fgzx_c(:,:,1) = i;
    i = fgzx_c(:,:,2); i(senszx) = par.senscolor(2);fgzx_c(:,:,2) = i;
    i = fgzx_c(:,:,3); i(senszx) = par.senscolor(3);fgzx_c(:,:,3) = i;
    fgazx(senszx) = 1;
end

if ~isempty(par.roi),
    for jj = 1:size(par.roi,4),
        r = par.roi(:,:,:,jj);

        sensxy = max(logical(r),[],3);
        sensxy = bwperim(sensxy); 
        i = fgxy_c(:,:,1); i(sensxy) = par.roicolor(1);fgxy_c(:,:,1) = i;
        i = fgxy_c(:,:,2); i(sensxy) = par.roicolor(2);fgxy_c(:,:,2) = i;
        i = fgxy_c(:,:,3); i(sensxy) = par.roicolor(3);fgxy_c(:,:,3) = i;
        fgaxy(sensxy) = 1;

        sensyz = squeeze(max(logical(r),[],1))';
        sensyz = bwperim(sensyz);
        i = fgyz_c(:,:,1); i(sensyz) = par.roicolor(1);fgyz_c(:,:,1) = i;
        i = fgyz_c(:,:,2); i(sensyz) = par.roicolor(2);fgyz_c(:,:,2) = i;
        i = fgyz_c(:,:,3); i(sensyz) = par.roicolor(3);fgyz_c(:,:,3) = i;
        fgayz(sensyz) = 1;

        senszx = squeeze(max(logical(r),[],2));
        senszx = bwperim(senszx);[x,y]=ind2sub(size(senszx), find(senszx));
        i = fgzx_c(:,:,1); i(senszx) = par.roicolor(1);fgzx_c(:,:,1) = i;
        i = fgzx_c(:,:,2); i(senszx) = par.roicolor(2);fgzx_c(:,:,2) = i;
        i = fgzx_c(:,:,3); i(senszx) = par.roicolor(3);fgzx_c(:,:,3) = i;
        fgazx(senszx) = 1;

    end
end

% Figure and Axis Positioning
rdim = [size(bgxy,1) size(bgxy,2) size(bgyz,1)]*par.scale;
f = figure('Color','black');
fpos = get(f,'Position');
ofpos = fpos;
fpos(3) = spac*3+rdim(2)+rdim(3);
fpos(4) = spac*3+rdim(1)+rdim(3);
fpos(1) = 10;%fpos(1)-(fpos(3)-ofpos(3));
fpos(2) = 10;%fpos(2)-(fpos(4)-ofpos(4));
set(f,'Position',fpos);
clear fpos ofpos;

ax1 = axes('Units','pixel','Position',[spac spac+spac+rdim(3) rdim(2) rdim(1)]);
ax2 = axes('Units','pixel','Position',[spac spac rdim(2) rdim(3)]);
ax3 = axes('Units','pixel','Position',[spac+spac+rdim(2) spac+spac+rdim(3) rdim(3) rdim(1)]);


axes(ax1); ih1=image(bgxy_c); axis off;  hold on; ih1=image(fgxy_c); set(ih1,'AlphaData',fgaxy);
axes(ax2); ih2=image(bgyz_c); axis off; hold on; ih2=image(fgyz_c); set(ih2,'AlphaData',fgayz);
axes(ax3); ih3=image(bgzx_c); axis off; hold on; ih3=image(fgzx_c); set(ih3,'AlphaData',fgazx);


% if ~isempty(par.fgimg2),
%     axes(ax1); ih1=image(fgxy2_c); set(ih1,'AlphaData',fgaxy2);
%     axes(ax2); ih2=image(fgyz2_c); set(ih2,'AlphaData',fgayz2);
%     axes(ax3); ih3=image(fgzx2_c); set(ih3,'AlphaData',fgazx2);
% end

% if ~isempty(par.sensvol),
%     sensxy = max(logical(par.sensvol),[],3);
%     sensxy = bwperim(sensxy);
%     axes(ax1);ih1 = image(ind2rgb(double(sensxy)*64,gray));set(ih1,'AlphaData',sensxy*par.tickalpha);
%     sensyz = squeeze(max(logical(par.sensvol),[],1))';
%     sensyz = bwperim(sensyz);
%     axes(ax2);ih2 = image(ind2rgb(double(sensyz)*64,gray));set(ih2,'AlphaData',sensyz*par.tickalpha);
%     senszx = squeeze(max(logical(par.sensvol),[],2));
%     senszx = bwperim(senszx);
%     axes(ax3);ih3 = image(ind2rgb(double(senszx)*64,gray));set(ih3,'AlphaData',senszx*par.tickalpha);
% end


% add axis lines
if prod(par.pos) > 0,
    axes(ax1);line([coord(1) coord(1)],[1 size(xy,1)],'Color','white');
    axes(ax1);line([1 size(xy,2)],[coord(2) coord(2)],'Color','white');
%     axes(ax2);line([coord(1) coord(1)],[1 size(yz,1)],'Color','white');
    axes(ax2);line([1 size(yz,2)],[coord(3) coord(3)],'Color','white');
%     axes(ax3);line([1 size(zx,2)],[coord(2) coord(2)],'Color','white');
    axes(ax3);line([coord(3) coord(3)],[1 size(zx,1)],'Color','white');
end


% Copy Frame
fr = getframe(f);
cimg = fr.cdata(:,:,:);
% cimg=hardcopy(gca,'-Dzbuffer', '-r0');
% cimg = cimg(end-10-n:end-10-1,13:n+12,:);
close(f);

