function [mask mask_soft] = cropusim2D(im, offset )
% Perform cropping of the ultrasound image


% 
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               OFFSET -  depth of the initial contour [expressed in % of total image height]
% 
%       OUTPUT DESCRIPTION:
%                IMNEW - output image .
% 
   sz = size(im);
 
  
M_x = [-1 0 1];              %% Mask gradient x direction
M_y = M_x';                  %% Mask gradient y direction
I = flipud(imresize(im,[512 512]));
G_x = convn(I,M_x,'same');                %% Gradient in x direction
G_y = convn(I,M_y,'same');                %% Gradient in y direction


% post-process Gradient in y direction
Base = G_y;
Base(1,:) = 0;
Base(end,:) = 0;
BW = Base>0.1*max(Base(:));
se = strel('disk',3);
BW = imerode(BW,se);
BW(1:512-round(offset*512),:) = 0;

% initial membrane line
for i = 1:512
   nonzero = median(find(BW(:,i)));
   y(i) = nonzero; 
end

if length(find(isnan(y))) > 0.3*length(y);
   mask = ones(size(im));
   mask_soft = mask;
   return;
end

y(find(isnan(y))) = mean(y(~isnan(y)));

% fitted membrane line 
x = [1:512];
coeffs = polyfit(x, y, 2);
% Get fitted values
fittedX = linspace(min(x), max(x), 512);
fittedY = polyval(coeffs, fittedX);

mask = zeros(512,512);
for i = 1:512
mask(1:round(fittedY(i)),i) = 1;
end

% soft mask
D = bwdist(mask);
D = mat2gray(D);
D = abs(D-1);
mask_soft = D;
mask_soft(find(mask_soft<1)) = exp(6.908.*mask_soft(find(mask_soft<1)))./1000;
mask_soft = mask_soft.^2;

% resize to the original size
mask = imresize(flipud(mask),sz);
mask_soft = imresize(flipud(mask_soft),sz);
   
   
end