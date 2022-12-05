function plotColorbar(pos,cmap,varargin)


par.cbmode = 1;             % alpha blending mode (if alpha map used)
par.cbarlabel = cell(1,3);
par.cbarbgcolor = [1 1 1];
par.cbarcolor = [0 0 0];
par.cbarfontsize = 8;
par.cbrotate = 0;
par.athres = [];

% Copy parameters from input struct
if numel(varargin) >= 1
    cpar = varargin{1};
    % transfer all fields to parameter array
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end
   
% pink colormap
cmap = getcmap(cmap);

% extend color map with alpha thresholds
if numel(par.athres) == 2,
    amap = zeros(1,64)';
    astart = round(par.athres(1)*64); if astart == 0, astart = 1; end;
    aend = round(par.athres(2)*64); if aend > 64, aend = 64; end;
    amap(astart:aend) = linspace(0,1,aend-astart+1);
    amap(aend:end) = 1;
   cmap = [cmap flipud(amap)];
end

if (ishandle(pos))
    axes(pos);
    set(pos,'Units','pixel');
    pos = uint16(get(pos,'Position'));
%     xp = pos;
%     pos(3) = 
%     set(gca,'Units','normalized');
else
    axes('Units','pixel','Position',pos);
end




%%
if (par.cbrotate)
    xvec = linspace(0,1,pos(3));
    yvec = ones(1,pos(4))';
else
    yvec = linspace(1,0,pos(4))';
    xvec = ones(1,pos(3));
end

if (size(cmap,2) == 4)
    cb = logical(checkerboard(4,ceil(pos(3)/8)+1,ceil(pos(4)/8)+1)); 
%     size(cb)
%     pos
    cb = cb(1:pos(3),1:pos(4))';
    image(ind2rgb(cb*64,gray));
    hold on;
end

% figure;
img_c=ind2rgb(round(yvec*xvec*64),cmap);
ih = image(img_c);
% axis image
axis off;

% if colorbar should get Alpha Blending
if (size(cmap,2) == 4)
    aimg = ones(pos(4),pos(3));     % start with opaque for all pixels
    % mode 1: show whole colorbar with alpha scaling
    if par.cbmode == 1
        aimg(:,:) = spline(linspace(0,1,64),cmap(:,4),linspace(0,1,pos(4)))'*xvec;
    % mode 2: show half of the colorbar with alpha scaling, half opaque
    elseif par.cbmode == 2        
        xvec = ones(1,ceil(pos(4)/2));
        aimg(floor(pos(4)/2)+1:end,:) = xvec'*spline(linspace(0,1,64),cmap(:,4),linspace(0,1,pos(3)));
    % mode 3: Smoothly blend alpha
    elseif par.cbmode == 3        
        avec = spline(linspace(0,1,64),cmap(:,4),linspace(0,1,pos(4)))';
        for j = 1:pos(4)
            aimg(j,:) = linspace(avec(j),1,pos(3));
        end
        if par.cbrotate
            aimg = aimg(end:-1:1,:);
        end
    end
    set(ih,'AlphaData',aimg);
end


%% labels
if par.cbarlabel{2} 
    h=text(double(pos(3))/2,5,par.cbarlabel{2},'rotation',90,'backgroundcolor',par.cbarbgcolor,'FontSize',par.cbarfontsize,'HorizontalAlign','right','Margin',0.5,'Color',par.cbarcolor);
end
if par.cbarlabel{1} 
    h=text(double(pos(3))/2,double(pos(4))-5,par.cbarlabel{1},'rotation',90,'backgroundcolor',par.cbarbgcolor,'FontSize',par.cbarfontsize,'HorizontalAlign','left','Margin',0.5,'Color',par.cbarcolor);
end
if numel(par.cbarlabel) > 2 && (numel(par.cbarlabel{3} > 0)),
    h=text(double(pos(3))/2,double(pos(4))/2,par.cbarlabel{3},'rotation',90,'backgroundcolor',par.cbarbgcolor,'FontSize',par.cbarfontsize,'HorizontalAlign','center','Margin',0.5,'Color',par.cbarcolor);
end