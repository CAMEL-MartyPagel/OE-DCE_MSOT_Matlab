function [amf,par] = msp_amfQLshr(Recon_orig,spec_orig,varargin)
% MSP_AMF  Adaptive Matched Filter ("AMF") unmixing with improved
% covariance matrix calculation
%   amf = msp_amfQLshr(Recon,spec)
%   amf = msp_amfQLshr(Recon,spec,par)

%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)

% Quasi-local covariance shrinkage: local covariance matrix is merged with a global cov. matrix
% Reference: Tzoumas et al. Statistical molecular target detection framework for multispectral optoacoustic tomography, IEEE TMI, 2016



par.thresh = 0;


% Copy parameters
if numel(varargin) >= 0,
    cpar = varargin{1};
    fx = fieldnames(cpar);
    for j = 1:numel(fx)
        par = setfield(par,fx{j},getfield(cpar,fx{j}));
    end
    clear cpar j fx;
end

[Recon_orig par wbar] = msp_preprocess(Recon_orig,par);

str = [];

% check if more than one spectrum is selected, then cancel
% TODO
if size(spec_orig,1) > 1
    errordlg('Please select only one spectrum for AMF');
    amf = [];
    return;
end
  
% % reduce to 2D
% svec = size(Recon);
% Recon = reshape(Recon,[prod(svec(1:end-1)) svec(end)]);

% check the wavelengths set
wavelengths = par.wavelengths;
req_wav = 700:10:900;
    flag_wav = 1;
    w_subset = intersect(wavelengths,req_wav);

    if length(w_subset)<length(wavelengths)        
        uiwait( warndlg('This algorithm was developed for wavelengths within [700-900] nm and multiple of 10. Some wavelengths of the current dataset will be discarded','modal'));
    end
    
    for i=1:length(w_subset) 
       [dummy,tmp] = min(abs(wavelengths - w_subset(i)));
       keep_idx(i) = tmp;
    end

    for i=1:length(w_subset) 
        ttmp = Recon_orig(:,:,:,:,1,keep_idx(i));
        Recon(:,:,:,:,1,i) = ttmp;
        spec(:,i) = spec_orig(:,keep_idx(i));
    end
 
   wavelengths = w_subset;
   %
   load('RSDF_data');  
   K = cov_gl;

%% AMF QL Shrickage

% Initialise Progress Bar
wbar = waitbar(0,'AMF QL Shr. Processing...' ) ;
count = 1;
for run=1:size(Recon,3)
    for slice=1:size(Recon,4)
        waitbar(count/(size(Recon,3)*size(Recon,4)),wbar);
      
        R = squeeze(Recon(:,:,run,slice,1,:));
        HM(:,:,1,:) = R;
        mixed = permute(HM,[4,3,1,2]);
        
        n1 = size(mixed,3);
        n2 = size(mixed,4);
        nm_sp = size(spec,1);
        A_est = spec';
        
        %Normalize Spectra
        for i=1:size(A_est,2)
            A_est(:,i) = A_est(:,i)./norm(A_est(:,i),2);
        end
        S = A_est;

        tmp = squeeze(mixed);
        x = tmp(:,:);
        [x, mixedmean] = remmean(x);
        psi1 = n1*n2; psi2 = 0;
         
        %% based on the current  wavelength sampling, recover the index of the prior
        if min(wavelengths)<700 || max(wavelengths)>900
            disp('No prior modeling available for current wavelengths');
            return;
        end
        w_prior = 700:10:900;
        for w_idx=1:length(wavelengths)
           [dummy idx_sampling(w_idx)] = min(abs(w_prior - wavelengths(w_idx))) ;
        end
     
        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G_init = G_init./norm(G_init,'fro');
        C_QL = C_QL./norm(C_QL,'fro');
               
        for kk=1:length(covAll)
            G_prior_i = covAll{kk}; G_prior_i_sampled = G_prior_i(idx_sampling,idx_sampling);
            G_prior_i_sampled = G_prior_i_sampled./norm(G_prior_i_sampled,'fro');
            G_dist_All(kk) = norm(G_init./norm(G_init,'fro')-G_prior_i_sampled,'fro');
        end
        
        a2 =  sort(G_dist_All);
        metric = a2(1);
        betha = 4*sqrt(metric); if betha>1 betha=1; end

        %Quasi local estimator of covariance matrix
        G_init = cov(x');
        [V D] = eig(K(idx_sampling,idx_sampling)); %K=V*D*V'
        C_QL = V*diag(diag(V'*G_init*V))*V';
        
        G = (1-betha)*G_init+betha*C_QL;
        G_inv = inv(G);    

        unmixed = (S'*G_inv*x);
        unmixed(unmixed<0) = 0;
        den = psi1*(S'*G_inv*S);
        unmixed = (unmixed.^2)./den;
        umx = reshape(unmixed,n1,n2);
        umx(umx<par.thresh) = 0;
%         umx = sqrt(umx);
        p = size(A_est,2);
        A = betha;
           
        amf(:,:,run,slice) = squeeze(umx);
%         if max(umx(:)>0)
%             umx = umx./max(umx(:));
%         end
     count = count + 1;     
    end
end
close(wbar);
   

