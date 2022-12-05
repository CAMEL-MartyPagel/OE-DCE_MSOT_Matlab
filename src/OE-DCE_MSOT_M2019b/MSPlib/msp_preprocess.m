function [Recon par wbar] = msp_preprocess(Recon,cpar)
% MSP_PREPROCESS  Common preprocessing for unmixing library
%   [Recon par] = msp_preprocess(Recon,par)
%
% Common preprocessing functions for MSP Processing:
% 
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - par: Parameters for Preprocessing
%   o imresize: Integer number for resolution reduction
%   o rejectneg: Reject Negatives from number of negative wls onwards (1)
%                (1 means if any wavelength is negative px is rejected,
%                 2 means rejection as soon as 2 or more wavelengths are -,
%                 0.33 means at least 33% of wls need to be <0 for reject)

% *** DEFAULTS ***
par.imresize = 0;
par.rejectneg = 1;
par.rejectlevel = 0;
par.is3D = 0;
par.progress = true;

% image sizes
svec = size(Recon);
nx = svec(1); ny = svec(2);
if (numel(svec) == 2), svec = [nx ny 1]; end;
par.initdim = svec;

% Copy parameters from input struct
% transfer all fields to parameter array
fx = fieldnames(cpar);
for j = 1:numel(fx)
    par = setfield(par,fx{j},getfield(cpar,fx{j}));
end
clear cpar j fx;

% Initialise waitbar
if par.progress,
    wbar = waitbar(0,'MSP Preprocessing...');
else 
    wbar = [];
end


%% *** Reject negatives
if par.rejectneg > 0,
    if par.rejectneg < 1, % given in percent
        par.rejectneg = ceil(par.rejectneg * svec(end));
    end
    if abs(par.rejectlevel) > 0 && abs(par.rejectlevel) < 1,
%         range = diff([min(Recon(:)) max(Recon(:))]);
        range = max(Recon(:));
        par.rejectlevel = par.rejectlevel * range;
    end
    
    Recon = reshape(Recon,[nx*ny*prod(svec(3:end-1)) svec(end)]);
    loc = Recon <= par.rejectlevel;
    negwl = uint8(sum(loc,2));
    reject = negwl >= par.rejectneg;
    loc = logical(reject*ones(1,svec(end)));
    par.rejectpx = nnz(reject);
    Recon(loc) = 0;
    
    Recon = reshape(Recon,[nx ny svec(3:end)]);    
end

%% *** IMAGE RESIZING
if (par.imresize > 0),
    if par.is3D,
        nz = svec(3);
        nimg = prod(svec(4:end));
        Recon = reshape(Recon,nx,ny,nz,nimg);

        nx = double(round(nx ./ par.imresize));
        ny = double(round(ny ./ par.imresize));
        nz = double(round(nz ./ par.imresize));
        ix = linspace(1,svec(1),nx);
        iy = linspace(1,svec(2),ny);
        iz = linspace(1,svec(3),nz);
        [gx gy gz] = meshgrid(ix,iy,iz);

        % get adjusted resolution
        if par.progress,wbar = waitbar(0,wbar,'Resizing Volumes...');end
        c = 0;
%         [nx ny nimg]
        R = zeros(nx,ny,nz,nimg);
        for sl = 1:nimg
           c = c+1;
           if par.progress,waitbar(c/nimg,wbar);end
           R(:,:,:,sl) = interp3(Recon(:,:,:,sl),gx,gy,gz);
        end
        Recon = reshape(R,[nx ny nz svec(4:end)]);
        clear R;
        
    % 2D case
    else
        nimg = prod(svec(3:end));
        Recon = reshape(Recon,nx,ny,nimg);

        nx = double(round(nx ./ par.imresize));
        ny = double(round(ny ./ par.imresize));

        % get adjusted resolution
        if par.progress,wbar = waitbar(0,wbar,'Resizing Images...');end;
        c = 0;
%         [nx ny nimg]
        R = zeros(nx,ny,nimg);
        for sl = 1:nimg
           c = c+1;
           if par.progress,waitbar(c/nimg,wbar);end
           R(:,:,sl) = imresize(Recon(:,:,sl),[nx ny]);
        end
        Recon = reshape(R,[nx ny svec(3:end)]);
        clear R;
    end
end     % imresize
close(wbar);
