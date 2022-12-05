function wl = WL_invert(wl,p_w,A_mat,n,t,n_proj,verbose)
% threshold to select interesting spots in the model relevant for
% reconstruction
over_fact=ceil(wl.over_fact*size(A_mat,1)/size(A_mat,2));

if (verbose) fprintf('Starting wl-domain matrix inversion...\n'); end
for kk=1:4^(wl.depth)
    if (verbose) fprintf('Wavelet packet %i of %i\n',kk,4^(wl.depth)); end;
    
    % this is the final wavelet inverse matrix R^i in eq. 24
    wl.Ri{kk} = calc_Ri(kk,n,n,wl.thres,wl.depth,wl.wl_name);

    % this is the wavelet transform of the model (eq. 24)
    % (the matrix D is the wavelet transform operator)
    Mi_w = transform_t2wl(A_mat*wl.Ri{kk},numel(t),n_proj,wl.depth,wl.wl_name);

    % find interesting places
    wl.places{kk} = choose_timeproj( Mi_w,over_fact,size(wl.Ri{kk},2)  ) ;

    % select places in both signal and image for eq. 26
    Mi_w_dash=Mi_w(wl.places{kk},:);
    clear Mi_w;

     % invert partial model matrix for this base
     [u,s,v]=svd(full(Mi_w_dash),'econ'); % economy version of SVD
     clear Mi_w_dash

     Mi_tmp = tsvd(u,s,v,wl.Sthres,max(s(:)));
     
     % threshold small elements in inverted model
     Mi_tmp (abs(Mi_tmp) < max(abs(Mi_tmp(:)))*wl.Athres) = 0;
     wl.Mi_w_dash_inv{kk}=sparse(Mi_tmp);
    clear s ss u v Mi_tmp;

end