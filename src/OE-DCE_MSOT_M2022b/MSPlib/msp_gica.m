function [umx par A W] = msp_gica(Recon,spec,varargin)
% MSP_GICA  Guided ICA unmixing
%   [gica par A W] = msp_gica(Recon,spec)
%   [gica par A W] = msp_gica(Recon,spec,par)
%
% Guided ICA unmixing using the FastICA library
% 
% Last dimension must be wavelengths, returns unmixed images in same
% dimensionality as the input recon
%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)
%   Specific to Guided ICA:
%   o mu: Step size (0.01)
%   o epsilon: Convergence (5e-5)
%   o verbose: Verbosity of algorithm (off)
%   o components: Number of components (number of input spectra)
%   o reapply: Do not optimize spec using ICA, but only multiply mean-free (0)

% check for ICA on path
if exist('fastica') ~= 2,
    error('FastICA not found in MATLAB path');
end;

par.components = size(spec,2);
par.mu = 0.01;
par.epsilon = 5e-5;
par.verbose = 'off';
par.method = 'gica';
par.reapply = 0;

% Copy parameters
if numel(varargin) > 0,
    cpar = varargin{1};
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

%% *** Reject negatives
% Only remember positions of negatives to not alter statistics, then zero
% after processing
svec = size(Recon);
nx = svec(1); ny = svec(2);
par.rejectlevel = 0;
if par.rejectneg > 0,
    if par.rejectneg < 1, % given in percent
        par.rejectneg = ceil(par.rejectneg * svec(end));
    end
    if abs(par.rejectlevel) > 0 && abs(par.rejectlevel) < 1,
%         range = diff([min(Recon(:)) max(Recon(:))]);
        range = max(Recon(:));
        par.rejectlevel = par.rejectlevel * range;
    end
    
    loc = sum(Recon < par.rejectlevel,6);       % count negative wavelengths per px
    loc = loc >= par.rejectneg;     % save location of negatives
   
end
par.rejectneg = 0;


[Recon par wbar] = msp_preprocess(Recon,par);


%% ICA
svec = size(Recon);
% Reshape to 2D vector
Recon = reshape(Recon,[prod(svec(1:end-1)) svec(end)]);

if size(spec,2) > size(Recon,2)
    umx = []; par = []; A = []; W = [];
    errordlg('Number of components > Number of wavelengths','MSP Error');
    return;
end

if (par.reapply == 0),
    info.factorICA = 1000000;
    mu = 0.01;
    eps = 5e-5;
    verb = 'on';
    [icasig, A, W] = fastica(Recon', info,'numOfIC', par.components , 'stabilization','on', 'initGuess', spec, 'mu', par.mu, 'approach', 'symm', 'epsilon', par.epsilon, 'verbose',par.verbose);
else
    par.components = size(spec,1);
    [Recon, mixedmean] = remmean(Recon');
    icasig = spec * Recon + (spec * mixedmean) * ones(1, prod(svec(1:end-1)));  
    W = spec;
    A = [];
end

icasig(:,loc) = 0;
% the unmixed array is converted back to a multi-D image
umx = reshape(icasig',[svec(1:end-1) par.components]);

%% Postprocessing
[umx par] = msp_postprocess(umx,par,wbar);

