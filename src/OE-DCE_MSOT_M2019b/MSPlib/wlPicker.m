function [ wl_optim, spec_optim, ind_optim ] = wlPicker( spec_full, wl_full, wl_fixed, N)

% this function automatically selects the wavelength set based on the criterion
% of the minimal singular value of the matrix of molar absorption at
% different wavelengths

% INPUT
% spec_full - initially given spectra of absorbers
% wl_full - initially given range of wavelengths
% wl_fixed - initially chosen wavelengths (fixed)
% N - number of desired wavelengths to be selected

% OUTPUT
% Spec_optim - molar absorption spectra at the selected wavelengths
% ind_optim - indices of the selected wavelengths from the full range

spec_optim = spec_full;
wl_optim = wl_full;
wlNum_desired = size(spec_optim,1);


while wlNum_desired > N %  until we reach the desired number of wavelengths
    [~, ind_fixed] = ismember(wl_fixed, wl_optim);
for i = 1:size(spec_optim,1)
spec_tmp = spec_optim;

if (ismember(i, ind_fixed) == 0);
spec_tmp(i,:) = []; % remove the each row iteratively

% compute singular values of the resulting matrix of absorption spectra
[U,S,V] = svd(spec_tmp);
if (size(S,1) == 1) || (size(S,2) == 1) 
    Smin(1,i) =  S(1,1);
else
Smin(1,i) =  min(diag(S)); % store the smallest singular value
end

end
end

[Smin_largest, ind_max] = max(Smin);
spec_optim(ind_max,:) = []; %  discard permanently the wavelengths removal of corresponding row of which led to the largest singular value
wl_optim(ind_max,:) = [];

wlNum_desired = size(spec_optim,1);
clear Smin ind_max Smin_largest
end


[~, ind_optim] = ismember(wl_optim, wl_full);
spec_optim = spec_full(ind_optim,:);


end

