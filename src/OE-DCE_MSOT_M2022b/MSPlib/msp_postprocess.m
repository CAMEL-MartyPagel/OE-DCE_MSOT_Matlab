function [Recon par] = msp_postprocess(Recon,cpar,wbar)
% MSP_POSTPROCESS  Common postprocessing for unmixing library
%   [Recon par] = msp_postprocess(Recon,par)
%
% Common postprocessing functions for MSP Processing:
% 
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - par: Parameters for Postprocessing
%   o imresize: Integer number for resolution reduction 

par = struct;
par.is3D = 0;

% image sizes
svec = size(Recon);
nx = svec(1); ny = svec(2);
if (numel(svec) == 2), svec = [nx ny 1]; end;


% Copy parameters from input struct
% transfer all fields to parameter array
fx = fieldnames(cpar);
for j = 1:numel(fx)
    par = setfield(par,fx{j},getfield(cpar,fx{j}));
end
clear cpar j fx;


% Initialise waitbar
if par.progress && ~isempty(wbar),
    wbar = waitbar(0,'MSP Postprocessing...' ) ;
end

%% *** IMAGE RESIZING
if (par.imresize > 0),
    if par.is3D,
        nz = svec(3);
        nimg = prod(svec(4:end));
        Recon = reshape(Recon,nx,ny,nz,nimg);

        % restore initial sizes
        nx = par.initdim(1);
        ny = par.initdim(2);
        nz = par.initdim(3);
        ix = linspace(1,svec(1),nx);
        iy = linspace(1,svec(2),ny);
        iz = linspace(1,svec(3),nz);
        [gx gy gz] = meshgrid(ix,iy,iz);

        % get adjusted resolution
        if par.progress,wbar = waitbar(0,wbar,'Restoring Volumes...');end
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

        % restore initial sizes
        nx = par.initdim(1);
        ny = par.initdim(2);

        % get adjusted resolution
        if par.progress,wbar = waitbar(0,wbar,'Restoring Image Resolution...');end
        c = 0;
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

%% close waitbar
if par.progress && ~isempty(wbar),
    close(wbar);
end