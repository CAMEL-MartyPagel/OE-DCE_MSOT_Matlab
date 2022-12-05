function filename = invertWL(A_mat,cpar)

% default settings
par = getWLdefaults;

% Input parameters overwrite defaults
% transfer all fields to parameter array
fx = fieldnames(cpar);
for j = 1:numel(fx)
    par = setfield(par,fx{j},getfield(cpar,fx{j}));
end
clear cpar j fx;

% Threshold of interesting places
over_fact=ceil(1.4*size(A_mat,1)/size(A_mat,2));

% Wavelet Filter Coefficients
[Lo_D,Hi_D] = wfilters(par.wl_name);


% invert using wavelet packets
% *************************************************************************
wbar = waitbar(0,'Estimated Time Left: Unknown','Name',...
    'Wavelet Domain Matrix Inversion','CreateCancelBtn',...
    'setappdata(gcbf,''canceling'',1)');

fprintf('Starting wl-domain matrix inversion...\n');
for kk=1:4^(par.depth_im)
    fprintf('Wavelet packet %i of %i\n',kk,4^(par.depth_im));

    if getappdata(wbar,'canceling')
        close(wbar);
        delete(wbar);
        error('User cancelled the inversion');
    end

%     profile clear;
%     profile -memory on;
    
    % this is the final R^i in eq. 2
    tic
    Ri{kk} = xvue_calc_Ri(kk,par.n,par.n,par.thres,par.depth_im,par.wl_name);
    % this is the wavelet transform of the model (eq. 24)
    % (the matrix D is the wavelet transform operator)
    Mi_w = mul_transform_t2wlcl(A_mat, Ri{kk},size(A_mat,1)/par.proj,par.proj,par.depth_proj,Lo_D, Hi_D);

    % find interesting places
    places{kk} = choose_timeproj( Mi_w,over_fact,size(Ri{kk},2)  ) ;

    % select places in both signal and image for eq. 26
    whos Mi_w A_mat
    Mi_w=Mi_w(places{kk},:);

    % invert partial model matrix for this base
    [u,s,v]=svd(Mi_w,'econ'); % economy version of SVD
    clear Mi_w

     % truncate SVD
     Mi_tmp = tsvd(u,s,v,par.Sthres,max(s(:)));

     % threshold small elements in inverted model
     Mi_tmp (abs(Mi_tmp) < max(abs(Mi_tmp(:)))*par.Athres) = 0;
     sparsity(kk) = nnz(Mi_tmp)./numel(Mi_tmp(:));
     Mi_w_dash_inv{kk}=sparse2(Mi_tmp);
     whos Mi_tmp Mi_w_dash_inv Ri places
     clear s ss u v Mi_tmp;
     
%      profreport;
%      profile off;

        % time estimation
     perpacket = toc;
     fprintf('  Execution time: %.1f s -- Sparsity of inv. model: %.3f\n',perpacket,1-sparsity(kk));
     est = round(perpacket*(4^(par.depth_im)-kk)); est_unit = 's';
     if (est > 120) est = round(est / 60); est_unit = 'min'; end;    
     wbar = waitbar(kk/(4^(par.depth_im)),wbar,sprintf( 'Estimated Time Left: %i%s', est, est_unit )) ;
     
end

filename = getWLfilename(par) ;

wbar = waitbar(1,wbar,sprintf( 'Writing to disk...' )) ;
saveMiRiPlaces(Mi_w_dash_inv, Ri, places, par.n, filename);

close(wbar);
delete(wbar);
