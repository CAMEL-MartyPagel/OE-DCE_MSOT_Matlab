function imnew = laplacian2D(im)
% Apply laplacian filter
%   DIFF_IM = ANISODIFF2D(IM, NUM_ITER, DELTA_T, KAPPA, OPTION) perfoms 
%   conventional anisotropic diffusion (Perona & Malik) upon a gray scale
%   image. A 2D network structure of 8 neighboring nodes is considered for 
%   diffusion conduction.
% 
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               
% 
%       OUTPUT DESCRIPTION:
%                IM - image .
% 
   M_lap = [0 1 0; 1 -4 1; 0 1 1];            %% Laplacian Mask
   imnew = convn(im, M_lap,'same');  
   imnew = medfilt2(imnew,[2 2]);

end