function imnew = blinddeconv2D(im, reszFctr, psf_sz, psf_iter, varargin )
% Perform blind deconvolution of the image with unknown PSF of the blur
% kernel

% 
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               psf_sz   - size of the PSF
%               psf_iter - number of iterations for PSF
% 
%       OUTPUT DESCRIPTION:
%                IMNEW - output image
% 
   
  if numel(varargin)>0
     psf = varargin{1};
  end

   sz = size(im);
   sznew = round(reszFctr*sz);
   y = imresize(im,sznew);
   y = mat2gray(y);

   
   
   sf = [psf_sz, psf_sz];       % size of the PSF
   maxiter = [psf_iter, 1];           % number of iterations for f and x
   n = 1;   
   
   % MAIN
   % intialize the true image x with the first y padded with a constant
   

    sy = size(y);
    sx = sy + sf - 1;
    sf2 = floor(sf/2);
    x = median(y(:)) * ones(sx);
    x(sf2(1)+(1:sy(1)), sf2(2)+(1:sy(2))) = y;

   h = waitbar(0,'Please wait...');
   
  for i = 1:n
      
    [x, f1] = mbd(x, y, maxiter, psf);
    waitbar(i / n);

%     %plot
%     figure(100); hold on;
%     subplot(131), imagesc((y)), title(sprintf('observed image y%d', i)); axis equal, axis tight; caxis([0 1]);
%     if exist('f', 'var')
%         subplot(132), imagesc(f), title(sprintf('estimated PSF f%d', i)); axis equal, axis tight; 
%     end
%     subplot(133), imagesc(x(sf2(1):end-sf2(1),sf2(1):end-sf2(1))), title(sprintf('estimated image x%d', i)); axis equal, axis tight; caxis([0 1]);
%     colormap gray
%     drawnow
%     pause
  
  end
  
    close(h) 
    imnew = (x(sf2(1):end-sf2(1),sf2(1):end-sf2(1)));
    imnew = imresize(imnew,sz);
    Imax = max(im(:));
    Imin = min(im(:));
    imnew = newrange(imnew,[Imin Imax]);
   
end