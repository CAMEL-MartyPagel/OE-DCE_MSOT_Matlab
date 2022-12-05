function [umx, info] = msp_pcaica(image_MOrg, info );
% =========================================================================
% [umx, info] = unmix_M_pcaica(image_MOrg, info);
% ------------------------------------------------------------------------
% unmixes a array of images with Independent Component Analysis (ICA) after
% having applied a Principal Component Analysis (PCA) to sort the
% components in the order of their importance
% ------------------------------------------------------------------------
% Input:
% ------
% image_MOrg.........| the matrix with the images (Npx, Npy, Npz, Nf, 
% ...................| Nw, Nrep, NNpos, Nrun)
% info...............| Information and error log struct
% .-cropmap .........| logical map (Nx,Ny) to crop the image '0' = crop, '1' = use
% .-pcaica.compNum...| number of desired output components (max: number of wl)
% .-pcaica.compValid.| number of valid components (max: compNum)
% ...................| returned by the ICA (default/init: number of wavelengths = compValid)
% .-spectra..........| returned output spectra. size Nw x compValid
% .-returnMsg........| exception handling information.
%
% Output:
% -------
% umx........................| an array of component (unmixed) images
% info.......................| definition see under input parameters
%
% Uses:
% -----
% fastica....................| Fast ICA calculation package
%
% ========================================================================
% DESCRIPTION
% the input is 8D array (Npx, Npy, Npz, Nf, Nw, Nrep, Npos, Nrun) where:
% Npx % the number of pixels in x axis
% Npy % the number of pixels in y axis
% Npz % the number of pixels in z axis
% Nf  % the number of frames
% Nw  % the number of wavelengths
% Nrep % the number of repetitions
% NNpos % the number of Npositions
% Nrun % the number of Nrun
% The ICA is performed simultaneously for all slices  Nf, Nrep, Npos and Nrun.
% More information on the Fast ICA algorithm can be found at:
% http://www.cis.hut.fi/projects/ica/fastica/
% ========================================================================
%

%% Get cropmap and #components from info struct
cropmap = logical(info.cropmap);
% Error handlin if # components is missing
if (~ isfield(info,'pcaica'))
    error('doMSOT:unmix_M_pcaica:infoStructIncomplete','Definition of info.pcaica is missing.');
end
if (isfield(info.pcaica,'compNum'))
    NoOfComps = info.pcaica.compNum;
else
    error('doMSOT:unmix_M_pcaica:infoStructIncomplete','Definition of info.pcaica.compNum is missing.');
end

%% some pre-processing
% Determine the image matrix dimensions
[ Npx, Npy, Npz, Nf, Nw, Nrep, Npos, Nrun]=size(image_MOrg); % read the dimensions

% Permute Matrix 
image_M =  permute( image_MOrg, [ 5 8 7 6 4 3 1 2 ] ) ;


%% perfor PCA and then ICA
% the image array is converted to a 2D array (Nw, rest_of_signal) to be
% used with the FastICa algorithm
image_M = image_M(:, :, :, :, :, :, cropmap);
image_M = reshape(image_M, Nw, Nrun*Npos*Nrep*Nf*Npz*sum(cropmap(:)));

% *** Discard negatives
% Only remember positions of negatives to not alter statistics, then zero
% after processing
if isfield(info,'discard_negative_values') && info.discard_negative_values,
    loc = sum(image_M < 0,1);       % count negative wavelengths per px
    loc = loc >= ceil(Nw*0.25);     % save location of negatives
else
    loc = false(size(image_M));
end



image_M2 = permute(image_M,[2 1]);
[A, pcasig, W] = princomp(image_M2);
pcasig =  permute(pcasig,[2 1]); % Factor 1000 removed, 13.03.2012, Till
%pcasig = 1000 * permute(pcasig,[2 1]);
%pcasig =  0.3 * permute(pcasig,[2 1]); ; %for debugging

[pcaicasig, A, W] = fastica(pcasig( 1:NoOfComps, :),info );
NRes=size(pcaicasig,1);

if (NRes < NoOfComps ) % give back empty frames and WL vectors
   pcaicasig = [pcaicasig; zeros((NoOfComps-NRes),size(pcaicasig,2))] ;
   A = [A zeros(size(A,1),(NoOfComps-NRes))] ;
   info.pcaica.compValid = NRes ;
   newWarning.identifier = 'doMSOT:unmix_M_ica:ICAcomponentsReduced';
   newWarning.message = ['fastica resulted in fewer output than input components. ',num2str(NRes),' instead of ',num2str(NoOfComps),' .Replace by zeros.'];
   newWarning.flag = 1;
   info.returnMsg = [info.returnMsg newWarning];
   if (info.verbose); warning(newWarning.identifier,['[doMSOT:unmix_M_ica] |',newWarning.message]); end
end

% make negatives zero
pcaicasig(:,loc) = 0;

umx=zeros(NoOfComps, Nrun, Npos, Nrep, Nf, Npz, Npx, Npy);
umx(:, :, :, :, :, :, cropmap) = reshape(pcaicasig,NoOfComps,Nrun, Npos, Nrep, Nf, Npz, sum(cropmap(:)));

% find spectra
% the ipA matrix does not contain the correct spectra. Theese are recovered
spectra = image_M * pinv( pcaicasig ) ;

% Commented for quantitative unmixing. 21.02.2012, Stefan & Till
% make spectra possitive
% corr = sign( median( spectra ) ) ;
% for ii = 1 : NoOfComps
%     umx( ii, :, :, :, :, :, :, : ) = umx( ii, :, :, :, :, :, :, : ) * corr(ii) ;
%     spectra( :, ii ) = spectra( :, ii ) * corr(ii) ;
%     spectra( :, ii ) = spectra( :, ii ) / max( spectra( :, ii ) ) ;
%     umx( ii, :, :, :, :, :, :, : ) = umx( ii, :, :, :, :, :, :, : ) / max( spectra( :, ii ) ) ;
% end


info.spectra = spectra;
% repermute matrix
umx = permute( umx, [7 8 6 5 1 4 3 2] ) ; 

