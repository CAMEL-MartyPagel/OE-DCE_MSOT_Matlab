function [psinv par] = msp_pinv(Recon,spec,varargin)
% MSP_PINV  Pseudo-Inverse ("Linear Regression") unmixing
%   psinv = msp_pinv(Recon,spec)
%   psinv = msp_pinv(Recon,spec,par)
%
% Pseudo-Inverse ("Linear Regression") unmixing using pinv(spec).
% 
% Last dimension must be wavelengths, returns unmixed images in same
% dimensionality in psinv
%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)

% preprocessing, if parameters are supplied
if numel(varargin) > 0,
    par = varargin{1};
    [Recon par wbar] = msp_preprocess(Recon,par);
else 
    par = struct;
end
par.method = 'pinv';

specP = pinv(spec');
svecX = size(Recon);
Recon = reshape(Recon,[prod(svecX(1:end-1)) svecX(end)]);
% psinv = zeros(size(Recon,1),size(specP,2));
% for j = 1:size(specP,2)
%     psinv(:,j) = Recon*specP(:,j);
% end
psinv = Recon*specP;
psinv = reshape(psinv,[svecX(1:end-1) size(specP,2)]);

if numel(varargin) > 0,
    [psinv par] = msp_postprocess(psinv,par,wbar);
end