function [img_c cax img] = imthres(img,varargin)
% IMTHRES
% Threshold images
%
% [img, cax] = IMTHRES(img)
%   Automatically threshold image with limits in cax
% [img, cax] = IMTHRES(img,cax)
%   Threshold image based on supplied limits
% [img_c, cax] = IMTHRES(img,cax,cmap)
%   threshold image based on supplied limits and applies colormap

% use supplied limits
cax = [];
if (numel(varargin) >= 1),
    cax = varargin{1};
end

cthres = [0.05 0.995];
if (numel(varargin) >= 3),
    cthres = varargin{3};
end


% automatic thresholding...
if nnz(cax > 0) < 2
    [ne cent] = hist(img(:),255);
    cum = ne*triu(ones(255,255),0); cum = cum ./ max(cum(:));
    crophigh = find(cum > cthres(2)); if isempty(crophigh), crophigh = 255; end 
    croplow = find(cum < cthres(1)); if isempty(croplow), croplow = 1; end
    if isempty(cax) || isnan(cax(1)), cax(1) = cent(croplow(end)); end
    if isempty(cax) || numel(cax) == 1 || isnan(cax(2)), cax(2) = cent(crophigh(1)); end
%     cax = [cent(croplow(end)) cent(crophigh(1))]; % use center as thres
end

% Colormap
cmap = [];
if (numel(varargin) >= 2 & ischar(varargin{2})),
    cmap = getcmap(varargin{2});
end

% threshold
img(img < cax(1)) = cax(1);
img(img > cax(2)) = cax(2);
img = img - cax(1);
img = img ./ (cax(2)-cax(1));

% filter with Frangi or LogGabor
if (numel(varargin) >= 4),
   filt = varargin{4};
   if strcmp(filt,'Frangi') || strcmp(filt,'LogGabor'),
      fprintf('Applying filter...\n');
   end
end
    
% exit if no colormap given
if (isempty(cmap))
    img_c = img;
    return;
end

% apply colormap
img_c = round(img * size(cmap,1));
img_c = ind2rgb(img_c,cmap);