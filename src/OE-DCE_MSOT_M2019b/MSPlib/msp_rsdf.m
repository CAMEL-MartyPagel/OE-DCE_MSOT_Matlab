function [rsdf,par] = msp_rsdf(Recon_orig,spec_orig,varargin)
% MSP_RSDF  Robust Statistical Detection Framework
%   rsdf = msp_rsdf(Recon,spec)
%   rsdf = msp_rsdf(Recon,spec,par)
%
%
%
% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)


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

%% RSDF
% Initialise Progress Bar
wbar = waitbar(0,'RSDF Processing...' ) ;
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
        a =  sort(G_dist_All);
        metric = a(1);
   
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
        tt = sort(unmixed);
        
        test = mean(tt(end-10:end));
        
        if test<4e-3 & betha<1
            v = find_t_DOF(x);
            unmixed = zeros(1,size(x,2));
            for j=1:size(x,2)
                unmixed(j) = sqrt((v-1)/((v-2)+x(:,j)'*G_inv*x(:,j)))*(S'*G_inv*x(:,j))/sqrt((S'*G_inv*S));
            end 
            unmixed(unmixed<0) = 0;
            umx = reshape(unmixed,n1,n2);
            umx = umx.^2;
            umx(umx<par.thresh) = 0;
            umx = umx.^2;
            umx = umx.^2;

        else
            unmixed = (S'*G_inv*x);
            unmixed(unmixed<0) = 0;
            den = psi1*(S'*G_inv*S);
            unmixed = (unmixed.^2)./den;
%             unmixed = sqrt(unmixed);
%             unmixed = sqrt(unmixed);
            umx = reshape(unmixed,n1,n2);
            umx(umx<par.thresh) = 0;
        end

        p = size(A_est,2);
        A = S;
           
        rsdf(:,:,run,slice) = squeeze(umx);
%         if max(umx(:)>0)
%             umx = umx./max(umx(:));
%         end
    count = count + 1;    
    end
end
close(wbar);

function v = find_t_DOF(x)
    d = Mahalanobis(x); 
    v_vec=3.5:0.5:10;
    cnt=1;
    eX = eprob(d); 
    c = eX.eprob(end:-1:1);
    srt_d = sort(d);
    %default -4
    Pi = logspace(0,-4,1000);
    for cnt=1:length(v_vec)
%         cnt
        v = v_vec(cnt);
        d_theor_t_correct = (v/((v-2)*size(x,1)))*d;
        exc_t = 1-fcdf(sort(d_theor_t_correct),size(x,1),v);
           
        exc_metric = 0;
        for i=1:length(Pi)
            [dummy idx1] = min(abs(exc_t - Pi(i)));
            [dummy idx2] = min(abs(c - Pi(i)));
            exc_metric = exc_metric + abs(srt_d(idx1) - srt_d(idx2));
        end
        exc_m(cnt) = exc_metric;
    end
    [dummy idx_min] = min(exc_m);
    v = v_vec(idx_min);

function [d] = Mahalanobis(mixed)

%format of mixed is lambda * x * y;   
% switch size(size(mixed))
%     case 3

S = size(mixed); 
M = S(1); N = S(2);
miu = zeros(M,1);
T = zeros(M,N);
G = zeros(M);

%remmean
for i = 1:M
    miu(i,1) = 1/N *sum(mixed(i,:));
    T(i,:) = mixed(i,:)-miu(i,1); 
end 

%corrected
for j = 1:N
    G = T(:,j) * T(:,j)' + G;
end 

G = G / (N-1);% figure,imagesc(G);
G_inv = inv(G);

d = zeros(N,1);
for j = 1:N
    d(j) = (mixed(:,j)- miu(:,1))'* G_inv * (mixed(:,j)- miu(:,1));
end


function [eX] = eprob(X)
% eprob: calculates the exceedance probability for n column vectors in the
%        array [m n] X, where m are the observations. The probability is 
%        output in percent. eX is output as a structure (see Output Arguments).
%
% Usage: eX = eprob(X);
%
% Input Arguments:
%
%   X - [m n] vector where m are the observations and n are the number of
%   datasets for which the exceedance probability is to be calculated. 
%   The size of m must be the same for all datasets.
%
% Output Arguments:
%
%   eX - structure array containing all output data
%   ex.data - input data X [m n]
%   ex.r - the number of rows, m
%   ex.c - the number of datasets (columns), n
%   ex.sort - X input data sorted in descending order
%   ex.rank - single column matrix of the sorted data rank
%   ex.eprob - calculated exceedance probability (rank/m+1)
%
% Example:
%
%   X = randn(1000,1) % create randomly distributed dataset with 1000 values
%   eX = eprob(X);
%
% Author: Jeff A. Tuhtan
% e-mail address: jtuhtan@gmail.com
% Release: 1.0
% Release date: 20.08.2010


Scap = 10; % active operational energy storage capacity
% Scap = StorCapPercent eX average annual generation
eX = struct;

eX.data = X;
eX.r = size(eX.data,1); % no. rows
eX.c = size(eX.data,2); % no. cols

eX.sort = sort(eX.data,'descend'); % sorts data in descending order
eX.rank = (1:eX.r)';
eX.eprob = zeros(eX.r,1);
eX.eprob = eX.rank./(eX.r+1);

function p = fcdf(x,v1,v2);
%FCDF   F cumulative distribution function.
%   P = FCDF(X,V1,V2) returns the F cumulative distribution function
%   with V1 and V2 degrees of freedom at the values in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

if nargin < 3, 
    error('Requires three input arguments.'); 
end

[errorcode x v1 v2] = distchck(3,x,v1,v2);

if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end

%   Initialize P to zero.
p = zeros(size(x));

k1 = find(v1 <= 0 | v2 <= 0 | round(v1) ~= v1 | round(v2) ~= v2);
tmp   = NaN;
p(k1) = tmp(ones(size(k1)));

% Compute P when X > 0.
k = find(x > 0 & ~(v1 <= 0 | v2 <= 0 | round(v1) ~= v1 | round(v2) ~= v2));
if any(k), 
% use A&S formula 26.6.2 to relate to incomplete beta function 
    xx = v2(k) ./ (v2(k) + v1(k) .* x(k));
    p(k) = 1 - betainc( xx, v2(k)/2, v1(k)/2 );
end

   

