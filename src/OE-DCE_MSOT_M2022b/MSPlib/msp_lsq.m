function [lsqinv R2 aR2 fi]= msp_lsq(Recon,spec,varargin)
% MSP_LSQ ("Linear Regression") unmixing using LSQ
%   [lsqinv R2 fi] = msp_lsq(Recon,spec)
%   [lsqinv R2 fi] = msp_lsq(Recon,spec,par)
%
% ("Linear Regression") unmixing using lsqnonneg.
% 
% Last dimension must be wavelengths, returns unmixed images in same
% dimensionality in lsqinv
%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)
% Return Values:
% - lsqinv: Matrix with unmixed values (concentrations)
% - R2: Matrtix with coefficients of determination (spectral confidence map)
% - aR2: Matrtix with adjusted coefficients of determination 
% - fi: Matrix with modeled values (fitted Recon values)

% preprocessing, if parameters are supplied
if numel(varargin) > 0,
    par = varargin{1};
    [Recon par wbar] = msp_preprocess(Recon,par);
else 
    par = struct;
end
par.method = 'lsq';

%store Recon size
svecX = size(Recon);

%reshape Recon
Recon = reshape(Recon,[prod(svecX(1:end-1)) svecX(end)]);

%perform linear regression
lsqinv = zeros([size(Recon,1),size(spec,2)]);
parfor i=1:size(Recon,1)
    
    [res,tmin]=lsqnonneg(spec,Recon(i,:)');
    lsqinv(i,:)=res;
end

%modeled values
fi=lsqinv*spec';

%mean of observed data 
ym=repmat(mean(Recon,2),1,size(spec,1));

%total sum of squares
SStot=sum((Recon-ym).^2,2);

%regression sum of squares
SSreg=sum((fi-ym).^2,2);

%residual sum of squares
SSres=sum((Recon-fi).^2,2);

%R squared
%R2=SSreg./SStot;
R2=1-SSres./SStot;

%adjusted R squared
n=size(spec,1); % # of data set values
p=size(spec,2);% # of regressors
aR2=1-(1-R2)*((n-1)/(n-p-1));

%reshape output
R2=reshape(R2,svecX(1:end-1));
aR2=reshape(aR2,svecX(1:end-1));
lsqinv = reshape(lsqinv,[svecX(1:end-1) size(spec,2)]);
fi = reshape(fi,svecX);

%postprocessing
if numel(varargin) > 0,
    [lsqinv par] = msp_postprocess(lsqinv,par,wbar);
    
end


