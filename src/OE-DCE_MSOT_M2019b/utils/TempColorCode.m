clc
%%
[fname, fpath] = uigetfile({'*.tif'; '*.*'});
info = imfinfo([fpath '\' fname]);
num_images = numel(info);
%%
clear B;
for k = 1:num_images;
    B(:,:,k) = imread([fpath '\' fname], k, 'Info', info);
end
[m, i] = max(B, [], 3);

Tmax = uint8(i ./ max(i(:)) .* 255);
Cmax = uint8(m ./ max(i(:)) .* 255);

% cd(fpath)

dot=strfind(fname,'.');
fname1 = fname(1:dot(end)-1);
imwrite(Tmax, sprintf('%s\\%s_Tmax.tif',fpath,fname1));
imwrite(Cmax, sprintf('%s\\%s_Cmax.tif',fpath,fname1));

% cd('C:\Users\alexander.lernbecher\Documents\MATLAB')