function [img, meta] = loadProdigyImage(filename)
%%%%% Binary image formation reading
%%%%% Created by Gency Jeng 07/23/2015
%%%%% All rights reserved by S-Sharp Corp. 
[fid,errmsg] = fopen(filename,'r');
if fid == -1,
    error('Error loading file %s: %s',filename,errmsg);
end

%%% Header
meta.frameIdx = fread(fid,1,'uint64');          % Frame index that are collected
meta.epochTime = fread(fid,1,'uint64');         % time stamp in format of Unix time
meta.timeStamp = epoch2datenum(meta.epochTime);  % timestamp in matlab format
meta.Width = fread(fid,1,'double');             % image width in mm
meta.Height = fread(fid,1,'double');            % image height in mm
meta.AxOffset = fread(fid,1,'double');          % Axial offset in mm
meta.dx = fread(fid,1,'double');                % beam interval in mm
meta.PixelWidth = fread(fid,1,'int');           % number of pixels in lateral direction 
meta.PixelHeight = fread(fid,1,'int');          % number of pixels in axial direction

%%% Image read
img = fread(fid,inf,'uint8');
img = reshape(img,meta.PixelHeight,meta.PixelWidth);

fclose(fid);
% figure
% imagesc(b)
% colormap(gray(256))
% title('Binary format')