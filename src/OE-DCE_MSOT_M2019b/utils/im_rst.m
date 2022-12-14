
function [rst]=im_rst(I1,x_shift,y_shift,angle)

% Function that performs Rotation, Scaling, Traslation of given image
% Input:
% I1: image being analysed
% x_shift: amount of shift along x-axis 
% y_shift: amount of shift along y-axis
% angle: angle of rotation counterclockwise [degrees]

I1_size = size(I1);
rst = zeros(I1_size);

scale = 1;        % scale factor
%angle = 0;         % rotation in degrees
%x_shift = 40;       % x translation (pixels)
%y_shift = -30;       % y translation (pixels)


% scale image
I_rst = imresize(I1,scale);
%rotate image
I_rst = imrotate(I_rst,angle,'bilinear');

size_larger = size(I_rst);

% figure(2),imshow(I_rst)

I_rst_shifted = zeros(size_larger);

if x_shift > 0
    x_start = x_shift+1;
    x_end = size_larger(2);
else
   x_start = 1;
   x_end = size_larger(2)+x_shift;
end

if y_shift > 0
    y_start = y_shift+1;
    y_end = size_larger(1);
else
   y_start = 1;
   y_end = size_larger(1)+y_shift;
end


for x = x_start:x_end
    for y=y_start:y_end
        I_rst_shifted(y,x) = I_rst(y-y_shift,x-x_shift);
    end
end


size_I = size(I_rst_shifted);
I_height = size_I(1);
I_width = size_I(2);

x_min = round(I_width/2 - I1_size(2)/2)+1;
y_min = round(I_height/2 - I1_size(1)/2)+1;

rst = imcrop(I_rst_shifted,[x_min y_min (I1_size(2)-1) (I1_size(1)-1)]);
rst = double(rst);

% imwrite(rst,'Im_rst_',num2str(scale),'_',num2str(angle),'_',num2str(x_shift),'_',num2str(y_shift),'bmp');