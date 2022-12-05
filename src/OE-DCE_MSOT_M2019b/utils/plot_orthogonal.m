function plot_orthogonal(varargin)

coord = [];
cax = [];
if (numel(size(varargin{1})) ~= 3)
    xy = varargin{1}; xy = xy(:,:,1);
    yz = varargin{2}; yz = yz(:,:,1)';
    zx = varargin{3}; zx = zx(:,:,1)';
else
    R = varargin{1};

    if (numel(varargin) > 1 && numel(varargin{2}) == 2)
        cax = varargin{2};
    end
        
    
    if (numel(varargin) > 2 && numel(varargin{3}) == 3)
        coord = varargin{3};
        xy = R(:,:,coord(3));
        zx = fliplr(squeeze(R(coord(1),:,:)));
        yz = fliplr(flipud(squeeze(R(:,coord(2),:))'));
    else
        xy = max(R,[],3);
        yz = squeeze(max(R,[],1))';
        zx = squeeze(max(R,[],2));
    end
end

spac = 10;
scale = 2;
rdim = [size(xy,1) size(xy,2) size(yz,1)]*scale;

f = figure;
fpos = get(f,'Position');
ofpos = fpos;
fpos(3) = spac*3+rdim(2)+rdim(3);
fpos(4) = spac*3+rdim(1)+rdim(3);
fpos(1) = fpos(1)-(fpos(3)-ofpos(3));
fpos(2) = fpos(2)-(fpos(4)-ofpos(4));
set(f,'Position',fpos);
clear fpos ofpos;

ax1 = axes('Units','pixel','Position',[spac spac+spac+rdim(3) rdim(2) rdim(1)]);
ax2 = axes('Units','pixel','Position',[spac spac rdim(2) rdim(3)]);
ax3 = axes('Units','pixel','Position',[spac+spac+rdim(2) spac+spac+rdim(3) rdim(3) rdim(1)]);


if isempty(cax)
    tmax = max([xy(:); yz(:); zx(:)]);
    tmin = min([xy(:); yz(:); zx(:)]);
    cax = [tmin tmax];
end

axes(ax1); imagesc(xy); caxis(cax); axis off; set(gca,'XColor','white');
axes(ax2); imagesc(flipud(yz)); caxis(cax); axis off;
axes(ax3); imagesc(fliplr(zx)); caxis(cax); axis off;

if ~isempty(coord)
    axes(ax1);line([coord(1) coord(1)],[1 size(xy,1)],'Color','white');
    axes(ax1);line([1 size(xy,2)],[coord(2) coord(2)],'Color','white');
%     axes(ax2);line([coord(1) coord(1)],[1 size(yz,1)],'Color','white');
    axes(ax2);line([1 size(yz,2)],[coord(3) coord(3)],'Color','white');
%     axes(ax3);line([1 size(zx,2)],[coord(2) coord(2)],'Color','white');
    axes(ax3);line([coord(3) coord(3)],[1 size(zx,1)],'Color','white');
end