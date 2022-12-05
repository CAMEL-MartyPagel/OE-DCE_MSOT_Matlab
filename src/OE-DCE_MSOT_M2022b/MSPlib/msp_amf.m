function [amf,par,bgmodel] = msp_amf(Recon,spec,varargin)
% MSP_AMF  Adaptive Matched Filter ("AMF") unmixing
%   amf = msp_pinv(Recon,spec)
%   amf = msp_pinv(Recon,spec,par)
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
%   o sphering: Data Sphering (1)
%   o suppbg: Background Suppression (0)


par.sphering = 1;
par.suppbg = 0;
par.bgmodel = [];

% Copy parameters
if numel(varargin) >= 0,
    cpar = varargin{1};
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

[Recon par wbar] = msp_preprocess(Recon,par);

str = [];

% check if more than one spectrum is selected, then cancel
% TODO
if size(spec,2) > 1
    errordlg('Please select only one spectrum for AMF');
    amf = [];
    return;
end
  
% reduce to 2D
svec = size(Recon);
Recon = reshape(Recon,[prod(svec(1:end-1)) svec(end)]);

% *** remove background
if par.suppbg
    bgi = mean(Recon');
    stdi = std(Recon');
    Recon(bgi < stdi,:) = 0;
end

x = Recon(:,:);
x = x';  % modified by EM 20170315 to correctly apply function remean() 


% Sphering of spectra and data
if par.sphering
    spec = spec ./ norm(spec);
    %str = [str sprintf('Removed mean from Spectra (%.1f)\n',Smean)];
%    Correct data sphering...
    [x, mixedmean] = remmean(x);
    str = [str sprintf('Removed mean from dataset (%.2f)\n',mean(mixedmean))];
end

% *** background modelling
if isempty(par.bgmodel),
    %G = cov(x);
    G = cov(x'); % modified by EM 20170315
    bgmodel = inv(G);
else
    bgmodel = par.bgmodel;
end
psi1 = svec(1)*svec(2);
psi2 = 0;  

% figure;
% subplot(1,3,1);
% imagesc(G);
% axis image;
% subplot(1,3,2);
% imagesc(G_inv);
% axis image;
% subplot(1,3,3);
% plot(diag(G_inv));

%Adaptive Matched Filter
%amf = (spec'*bgmodel*x');
amf = (spec'*bgmodel*x); % modified by EM 20170315
% figure;imagesc(reshape(unmixed,svec(1:2)))
amf(find(amf<0)) = 0;
den = spec'*bgmodel*spec;
str = [str sprintf('Denominator: %.2e\n',den)];
amf = (amf.^2)/den;
% unmixed = (unmixed.^2);
%
%CFAR property: Global threshold around 4e-4!
%3.8652e-004
thresh = 3.8652e-004;
%Calculate Pfa from Kelly's formula:
% Pfa = (1-thresh).^(psi1+1-size(x,2))*1e18
%         Wrong_pixels = Pfa*40000./100
%        unmixed(find(unmixed<thresh)) = 0;
%                 unmixed(find(unmixed>thresh)) = 1;

amf = reshape(amf,svec(1:end-1)); 
[amf par] = msp_postprocess(amf,par,wbar);
