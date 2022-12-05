function [ax1 ax2] = mip_2D(umx,varargin)
% MIP_2D
%   Create and Plot MIP of multi-slice dataset
%
%   [ax1 ax2] = MIP_2D(umx[,res,cmap])
% Arguments:
% - umx: 3D image matrix (x y z)
% - res: Resolution (aspect ratio)

ax1 = []; ax2 = [];

%% *** INPUT PARAMETERS
% Image Matrix
umx = squeeze(umx);
svec = size(umx);
if (numel(svec) > 3),
    error('Input has more than three dimensions');
end;
% Resolution
res = [ 1, 1 ];
if (numel(varargin) >= 1)
    res = varargin{1};
    if (numel(res) == 1), res = [res res]; end
end
if (res(1) ~= 1), res = res ./ res(1); end
% Colormap
cmap = gray(255);
if (numel(varargin) >= 2)
    cmap = varargin{2};
end




%% *** CREATE MIPS
mip_xy = max(umx,[],3);
mip_z = squeeze(max(umx,[],2));
[ny nz] = size(mip_z);
[Y,Z] = meshgrid(linspace(1,nz,round(nz*res(2))),linspace(1,ny,ny));
mip_z = interp2(double(mip_z),Y,Z);

% val = [mip_xy(:); mip_z(:)];
val = [mip_xy(:)];      % only threshold based on xy mip
[ne cent] = hist(val,100);
cum = ne*triu(ones(100,100),0); cum = cum ./ max(cum(:));
crophigh = find(cum > 0.995); % find elements exceeding % threshold
croplow = find(cum < 0.01); if isempty(croplow), croplow = 1; end
cax = [cent(croplow(end)) cent(crophigh(1))]; % use center as thres





%% *** CREATE PLOTS
f = figure;
whitebg(f);

ax1 = axes('Units','pixel','Position',[5 5 size(mip_xy,2) size(mip_xy,1)]);
image(imthres(mip_xy,cax,cmap));
axis off;

ax2 = axes('Units','pixel','Position',[5+5+size(mip_xy,2) 5 size(mip_z,2) size(mip_z,1)]);
image(imthres(mip_z,cax,cmap));
axis off;
