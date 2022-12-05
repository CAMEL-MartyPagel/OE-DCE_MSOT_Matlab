function [SA_NonLin, outer_mask] = msp_emsot_sO2(datapath, wavelengths, Recon, Mask, image_width, partMask)

close all;
clear U_d_container_norm;
clc

addpath('P:\_Steffen\software\eMSOT\')
addpath('P:\_Steffen\software\eMSOT\LightPropagation\2DFiniteVolumeMethod\Fluence_Simulations\fluence\fvm')

%load([datapath 'wavelengths']);

spectra = LoadSpectra(datapath,wavelengths);

%% Generate the eigenspectra model
method = 'analytical';
dim = 3;
[M, Model, ~, ~] = eigenspectra(datapath,method,dim,wavelengths);

Eig1 = Model(:,1)';
Eig2 = Model(:,2)';
Eig3 = Model(:,3)';

%% Compute prior values for m1, m2 and m3 based on a 3D or 2D light fluence model and first linear umx results
obj_mask = Mask;
MSPimage = squeeze(Recon(:,:,:));
SNR_thresh = 0;
[prior, ~] = createPriors2D(MSPimage,datapath,wavelengths,Mask,image_width,SNR_thresh);
%     [prior SA_prior_c] = createPriors2D_Adv(MSPimage,datapath,wavelengths,Mask,image_width,SNR_thresh);


%% STEP 6: Semi-Automatic grid creation on a selected tissue area
% An area on the image is segmented for MSOT application. The area must be
% of high intensity (high SNR) in a well reconstructed part of the image.
% The area must begin from tissue surface and go deep.
close all;

Recon_MB = Recon;


% area selection 

% % % % % radius_hole = 0;%ceil(min(size(obj_mask))/30);
% % % % % center_hole_W = size(obj_mask,1)/2;
% % % % % center_hole_H = size(obj_mask,2)/2;
% % % % % [W,H] = meshgrid(1:size(obj_mask,1),1:size(obj_mask,2));
% % % % % center_hole_mask = sqrt((W-center_hole_W).^2 + (H-center_hole_H).^2) > radius_hole;

outer_mask = obj_mask & partMask; % & center_hole_mask;
num_lines = 1 + ceil(min(size(obj_mask))/20); % minimum 2
num_points = 1 + ceil(min(size(obj_mask))/20); % minimum 2
[ X_n, Y_n, line_points_n, weight_rev, ~, ~] = grid_auto( Mask, outer_mask, num_lines, num_points , M, Model, spectra, Recon_MB);

% X_n, Y_n are the coordinates of the grid points. Each line of X_n/Y_n
% corresponds to a radial line of the grid. weight_rev contains 0 or 1
% depending on whether the grid point will be considered in the inversion

%% STEP 7: Prepare for Inversion
% At this point we should have: Recon, Recon_MB, Mask_3D, Mask, prior,
% param_save_path, phantom_name,datapath, Nx, Ny, Ns, Nw, wavelengths, M,
% Model, Eig1, Eig2, Eig3, outer_mask, num_lines, num_points, X_n, Y_n, line_points_n, weight_rev, c_x, c_y

close all;
clear points;

% points: num of lines x num of points x (X,Y)
% access in the image domain is image(Y,X).
points = line_points_n;

clear measured_s;
for l=1:size(line_points_n,1)
    for p=1:size(line_points_n,2)
        % measured_s: optoacoustic spectra corresponding to grid points
        measured_s(l,p,:) = squeeze(Recon_MB(points(l,p,2),points(l,p,1),:))./norm(squeeze(Recon_MB(points(l,p,2),points(l,p,1),:)));
        % priors: m1 and m3 prior values corresponding to the grid points.
        % The inverse order (p,l) instead of (l,p) is selected here for
        % avoiding extra permutations.
        prior_m1_mean(p,l) = prior(points(l,p,2),points(l,p,1),1);
        prior_m3_mean(p,l) = prior(points(l,p,2),points(l,p,1),3);
    end
end


%Check for NaN points in the measured_s and update the weights with zero
%values in all such cases.
tmp = isnan(measured_s);
if sum(tmp(:)) > 0
   warning('Data contain NaN values due to division with zero norm. Measured spectra and weights are accordingly updated') 
end
nan_idx = zeros(size(tmp,1),size(tmp,2));
measured_s(tmp) = 0;
for ii=1:size(tmp,3)
    nan_idx = nan_idx | tmp(:,:,ii);
end
weight_rev(nan_idx) = 0;

%% STEP 8: Select inversion parameters and apply Inversion

% Compute the maximum depth ( in cm) of each radial line of the grid for 
% enforcing depth dependent constraints (constrains of search space).
dist_depth = zeros(1,num_lines);
for i=1:num_lines
   dist_depth(i) = sqrt((X_n(i,1) - X_n(i,end)).^2 + (Y_n(i,1) - Y_n(i,end)).^2)*(100*image_width/size(Recon,1)); 
end

% Upper and lower limits of search space (with respect to the prior m1, m3) 
% identified ad hoc from the comparison of the prior to the actual values 
% in simulations of uniform optical properties and sO2
prior_m1_std_up = zeros(num_points,num_points);
for l=1:size(line_points_n,1)
   prior_m1_std_up(:,l) = linspace(0.02, 0.25*(dist_depth(l)/0.9), size(line_points_n,2)); 
end
prior_m1_std_low = zeros(num_points,num_points);
for l=1:size(line_points_n,1)
   prior_m1_std_low(:,l) = linspace(0.02, 0.18*(dist_depth(l)/0.9), size(line_points_n,2)); 
end
prior_m3_std_up = zeros(num_points,num_points);
for l=1:size(line_points_n,1)
   prior_m3_std_up(:,l) = linspace(0, 0.06*(dist_depth(l)/0.9), size(line_points_n,2)); 
end
prior_m3_std_low = zeros(num_points,num_points);
for l=1:size(line_points_n,1)
   prior_m3_std_low(:,l) = linspace(0, 0.06*(dist_depth(l)/0.9), size(line_points_n,2)); 
end


coords = line_points_n;                             % coordinates of the grid points
method = 'LocalOptimizerVectorNewInversionMAP';     % name of the optimizer
lambda = 0.3;                                       % value of Lagrange parameters for regularization

use_linear = 'on';                                  % initialization of sO2 with the linear unmixing result
constrain = 0;                                      % contrain to only m1<0 (contrain=1) or m1>0 (constrain=-1). This is not used int he current version   

% Update the global optimization limits based on the grid depth
lims = Compute_lims(wavelengths,datapath, mean(dist_depth));

% Apply Inversion
tic;                
[~, optans] = InversionArcRevisionMAP5(M, Model, lims, spectra, measured_s, coords, method, ...
lambda, constrain, use_linear, prior_m1_mean(:)', prior_m1_std_up(:)', prior_m1_std_low(:)', prior_m3_mean(:)', prior_m3_std_up(:)', prior_m3_std_low(:)', weight_rev,image_width,size(Recon,1));  
toc;

%% STEP 9: Interpolate fluence values for intermediate grid points, correct for fluence and compute sO2.

[Recon_MB_corrected, ~, ~, ~, ~, ~]  = interpolate_Fluence(X_n, Y_n, num_points, num_lines,optans, Recon_MB, M, Eig1, Eig2, Eig3, outer_mask);                
[SA_NonLin, error_NL, ~] = calc_sat(Recon_MB_corrected, spectra, wavelengths, 0); 
error_NL = sqrt(error_NL);

% Pixels with a fitting residual larger than 0.15 are not considered
% reliable
outer_mask(error_NL>0.15) = 0;

% The extreme values of -1 and 1 define the colormap visualization
SA_NonLin(~outer_mask) = -1; %SA_NonLin(end,end) = 1; 

