function [ wl_optim, spec_optim, ind_optim, cofa_mat_out ] = wlPicker2( spec_full, wl_full, wl_fixed, N, crit)
% WLPICKER2 Wavelength picker based on minimizing the angle between spectral vectors  
%
% INPUT
% spec_full - initially given spectra of absorbers
% wl_full - initially given range of wavelengths
% wl_fixed - initially chosen wavelengths (fixed)
% N - number of additional wavelengths to be picked
% crit - criterion to chose wls.
%        Use 'dist' to find a small angle, while having the wls equally
%        distributed, and covering the whole wl range.
%        Use 'minangle' to search for minimal angle between spectral vectors.
%    
% OUTPUT
% wl_optim - selected wavelengths from the full range
% spec_optim - spectra at the selected wavelengths
% ind_optim - indices of the selected wavelengths from the full range
% cofa_mat_out - matrix with cosine of the angles between the selected spectral vectors
%
% version 0.4 AU151117

% calculate number of possible WL combinations
comb_num=nchoosek(length(wl_full)-length(wl_fixed),N);

% wl index
wl_ind=[1:(length(wl_full))];
wl_ind_fixed=[];
for i=1:length(wl_fixed)
    wl_ind_fixed=[wl_ind_fixed,find(wl_full==wl_fixed(i))];
    wl_ind(find(wl_full==wl_fixed(i)))=[];
end

if comb_num<10^7 % check if number of combinations is below threshold
    comb0=nchoosek(wl_ind,N); %calculate all combinations
    comb=sort([comb0,repmat(wl_ind_fixed,size(comb0,1),1)],2); %add fixed wavelengths and sort
else
    'Too many combinations'
end

% intit variables
cofa_mat=zeros(size(spec_full,2),size(spec_full,2),comb_num);
cofa_norm0=zeros(comb_num,1);
equidist_measure=zeros(comb_num,1);
wl_sel_std=zeros(comb_num,1);

%loop over all combinations
parfor i=1:comb_num
       
    % spectra with selected wavelength
    spec_temp=spec_full(comb(i,:),:);
    
    % normalize
    for j=1:size(spec_temp,2)
        spec_temp(:,j)=spec_temp(:,j)/norm(spec_temp(:,j));
    end
    
    % cosine of angles (cofa) between spectral vectors
    cofa=spec_temp'*spec_temp;

    % remove diagonal elements and select upper triagonal values
    cofa=triu(cofa.*(~diag(ones(size(cofa,1),1)))); 
    cofa_mat(:,:,i)=cofa;

    % norm
    cofa_norm0(i)=norm(cofa(:));  
    
    % wl distance/distribution
    wl_sel=wl_full(comb(i,:));
    dist_mat=repmat(wl_sel,1,length(wl_sel))- repmat(wl_sel,1,length(wl_sel))'; % distances between wl
    equidist_measure(i)=std(dist_mat(dist_mat>0)); %std of distances
    
    % spread
    wl_sel_std(i)=std(wl_sel);
   
end

% parameter to be minimal
cofa_norm=cofa_norm0/max(cofa_norm0(:));
equidist_measure=equidist_measure/max(equidist_measure(:));
wl_sel_std=wl_sel_std/max(wl_sel_std(:));

if strcmp(crit,'dist')==1
    %finding a small angle, while having the wls equally distributed, and
    %covering the whole wl range
    min_par=cofa_norm.*equidist_measure./wl_sel_std; 
elseif strcmp(crit,'minangle')==1
    min_par=cofa_norm; %finding minimal angle
end

% index of combination of choice
comb_ind=find(min_par==min(min_par));

%ouput
ind_optim=comb(comb_ind,:);
wl_optim=wl_full(ind_optim);
spec_optim=spec_full(ind_optim,:);
cofa_mat_out=(cofa_mat(:,:,comb_ind));
 
end
