function [CosTheta fi] = msp_calcCos(Recon,U,spec,varargin)
% MSP_CALCCOS  Compute coefficients of determination (spectral confidence map)

% Parameters:
% - Recon: Reconstructed images, 3-XD matrix with wavelength as last dim
% - U: Matrix with unmixed values (concentrations)
% - spec: Spectra to unmix for
% - par: (optional) Parameters for Unmixing and Preprocessing (see
%                   msp_preprocess)
% Return Values:
% - R2: Matrtix with coefficients of determination (spectral confidence map)
% - aR2: Matrtix with adjusted coefficients of determination 
% - fi: Matrix with modeled values (fitted Recon values)
 
%store Recon size
svecX = size(Recon);
Recon0 = Recon;

%reshape Recon
Recon = reshape(Recon,[prod(svecX(1:end-1)) svecX(end)]);

%demean
ReconM = mean(Recon,2);
Recon_demean = Recon-repmat(ReconM,1,size(spec,1));

Usz = size(U);

if size(spec,2) == 1
    U = reshape(U,prod(Usz), 1 );
else
    U = reshape(U,prod(Usz(1:end-1)), Usz(end) );
end

%modeled values
fi=U*spec';
%demean
fiM = mean(fi,2);
fi_demean = fi-repmat(fiM,1,size(spec,1));
factor = ReconM./fiM;
factor(factor==-Inf) = 1;
factor(factor==Inf) = 1;
fi_demean1 = fi_demean.*repmat(factor,1,size(spec,1));

%% Cosine
CosTheta = (dot(Recon_demean',fi_demean1'))'./(vecnorm(Recon_demean,2).*vecnorm(fi_demean1,2));

%reshape output
CosTheta=reshape(CosTheta,svecX(1:end-1));
fi = reshape(fi,svecX);

% %%
% fi2 = reshape(fi,svecX );
% smod = squeeze(fi2(100,100,1,1,1,:));
% sreal = squeeze(Recon0(100,100,1,1,1,:));
% shb = spec(:,1);
% shbo2 = spec(:,2);
% srealM = mean(sreal);
% smodM = mean(smod);
% sreal_demean = sreal-srealM;
% smod_demean = smod-smodM;
% 
% factor = srealM/smodM;
% smod_demean1 = smod_demean*factor;
% wl = [680:5:715,730,760,800,850];
% 
% figure; plot(sreal); hold on; plot(smod)
% axis tight
% grid on
% axis square
% legend('Measured spectrum','Fitted spectrum');
% 
% figure; plot(sreal_demean); hold on; plot(smod_demean1)
% axis tight
% grid on
% axis square
% legend('Measured spectrum (de-meaned)','Fitted spectrum (de-meaned)');
% 
% SStot=sum((sreal'-srealM').^2,2);
% SSres=sum((sreal'-smod').^2,2);
% R2=1-SSres./SStot
% 
% SStot=sum((sreal'-srealM').^2,2);
% SSres=sum((sreal_demean'-smod_demean1').^2,2);
% R2=1-SSres./SStot
% 
% CosTheta_sreal_smod = dot(sreal_demean,smod_demean1)/(norm(sreal_demean)*norm(smod_demean1));
% 
% shbM = mean(shb);
% shbo2M = mean(shbo2);
% shb_demean = shb-shbM;
% shbo2_demean = shbo2-shbo2M;
% 
% factor = shbM/smodM;
% smod_demean1 = smod_demean*factor;
% CosTheta_shb_smod = dot(shb_demean,smod_demean1)/(norm(shb_demean)*norm(smod_demean1));
% figure; plot(smod_demean1)
% hold on; plot(shb_demean)
% axis tight
% grid on
% axis square
% legend('Fitted spectrum (de-meaned)','Hb spectrum (de-meaned)');
% 
% factor = shbo2M/smodM;
% smod_demean1 = smod_demean*factor;
% CosTheta_shbo2_smod = dot(shbo2_demean,smod_demean1)/(norm(shbo2_demean)*norm(smod_demean1));
% figure; plot(smod_demean1)
% hold on; plot(shbo2_demean)
% axis tight
% grid on
% axis square
% legend('Fitted spectrum (de-meaned)','HbO2 spectrum (de-meaned)');


% 
% %adjusted R squared
% n=size(spec,1); % # of data set values
% p=size(spec,2);% # of regressors
% aR2=1-(1-R2)*((n-1)/(n-p-1));
% 
% %reshape output
% R2=reshape(R2,svecX(1:end-1));
% aR2=reshape(aR2,svecX(1:end-1));
% fi = reshape(fi,svecX);
    
end
