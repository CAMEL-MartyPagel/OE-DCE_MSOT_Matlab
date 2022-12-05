function imnew = filterUSImage(im, method, par, varargin )
% Apply different image enhancement technique to US image


% 
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               METHOD   - method name (char) 
%               PAR - settings structure
% 
%       OUTPUT DESCRIPTION:
%               IMNEW - output image 
% 
   
mask = par.mask;
 
im(isnan(im)) = 0;
imnew = im;

if strcmp(method,'Anisotropic filter')
    
    % PARAMETERS
   num_iter = par.num_iter;
   kappa = par.kappa;
   
   % MAIN PROCESSING STEP
   imnew = anisodiff2D(im, num_iter, 1/7, kappa, 1);

elseif strcmp(method,'Laplacian filter')
   
   % MAIN PROCESSING STEP
   imnew = laplacian2D(im);
   
elseif strcmp(method,'Sigmoid mapping')  
   % PARAMETERS
   a = par.a;
   b = par.b;
   Imax = max(im(:));
   Imin = min(im(:));
   
   % MAIN PROCESSING STEP
   imnew = sigmoid2D(im, a, b);
   
elseif strcmp(method,'Deconvolution') 
    
   % PARAMETERS
   reszFctr = 0.8;
   psf_sz = par.psf_sz;
   psf_iter = par.psf_iter;

   if psf_sz == 50
    load('psf.mat');
   else
    psf = [];
   end
 
   % MAIN PROCESSING STEP
   imnew = blinddeconv2D(im, reszFctr, psf_sz, psf_iter, psf );
 
elseif strcmp(method,'Cropping') 
  imnew = im.*mask;

end

if ~strcmp(method,'Cropping') && ~isempty(mask)
    imnew = imnew.*mask;
end

   
end