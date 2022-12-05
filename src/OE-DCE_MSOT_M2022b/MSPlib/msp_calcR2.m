function [R2 aR2 dR2 cos fi] = msp_calcR2(Recon,U,spec,varargin)
% MSP_CALCR2  Compute coefficients of determination (spectral confidence map)

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

% normalize & demean
ReconM = mean(Recon,2);
Recon_demean = Recon-repmat(ReconM,1,size(spec,1));
Recon_norm = Recon ./norm(Recon,2);
Recon_demean1 = remmean(Recon_norm);
%Recon_demean1 = Recon_demean./repmat(max(Recon_demean,[],2),1,size(spec,1)); % normalize

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
fi_norm = fi ./ norm(fi,2);
fi_demean1 = remmean(fi_norm);
%fi_demean1 = fi_demean./repmat(max(fi_demean,[],2),1,size(spec,1)); % normalize

for i=1:size(spec,2)
    newSpec = spec(:,i) ./ norm(spec(:,i));
    [newSpec Smean] = remmean(newSpec);
    spec(:,i) = newSpec;
end

% %% spectral analysis for debugging only
% fi_tmp = reshape(fi,[332 332 25]);
% Recon_demean1_tmp = reshape(Recon_demean1,[332 332 25]);
% fi_demean1_tmp = reshape(fi_demean1,[332 332 25]);
% 
% wl = [660:10:900];
% figure; imagesc(Recon0(:,:,1,1,1,15));
% axis square;  axis off; colormap gray;
% 
% figure; plot(wl,squeeze(Recon0(140,300,1,1,1,:))); grid on; axis square; axis tight
% hold on; plot(wl,squeeze(fi_tmp(140,300,:))); grid on; axis square; axis tight
% legend('Measured','Fitted');
% R2(140,300)
% 
% figure; plot(wl,squeeze(Recon_demean1_tmp(140,300,:))); grid on; axis square; axis tight
% hold on; plot(wl,squeeze(fi_demean1_tmp(140,300,:))); grid on; axis square; axis tight
% legend('Measured','Fitted');
% dR2(140,300)
% 
% figure; plot(wl,squeeze(Recon0(155,240,1,1,1,:))); grid on; axis square; axis tight
% hold on; plot(wl,squeeze(fi_tmp(155,240,:))); grid on; axis square; axis tight
% legend('Measured','Fitted');
% R2(155,240)
% 
% figure; plot(wl,squeeze(Recon_demean1_tmp(155,240,:))); grid on; axis square; axis tight
% hold on; plot(wl,squeeze(fi_demean1_tmp(155,240,:))); grid on; axis square; axis tight
% legend('Measured','Fitted');
% dR2(155,240)

%% R2 V2

%total sum of squares
SStot=sum((Recon_demean1).^2,2);

%residual sum of squares
SSres=sum((Recon_demean1-fi_demean1).^2,2);

% R squared
dR2=1-(SSres./SStot);

%reshape output
dR2=reshape(dR2,svecX(1:end-1));

% % re-scale
% dR2Max = repmat(max(max(dR2,[],1),[],2),332,332);
% dR2 = dR2./dR2Max;

%% Cosine
cos = (dot(Recon_demean1',fi_demean1'))'./(vecnorm(Recon_demean1,2).*vecnorm(fi_demean1,2));

%reshape output
cos=reshape(cos,svecX(1:end-1));

% % re-scale
% cosMax = repmat(max(max(cos,[],1),[],2),332,332);
% cos = cos./cosMax;

%%

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
fi = reshape(fi,svecX);
    
end
