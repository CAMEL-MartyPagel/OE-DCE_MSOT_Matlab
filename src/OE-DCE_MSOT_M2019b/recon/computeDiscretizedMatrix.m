function A_mat = computeDiscretizedMatrix( c0, n, image_width, t, r_sensor, angle_sensor, n_angles )

nn = n * n ; % number of columns of the matrix
lt = length( t ) ;
n_rows = lt * length( angle_sensor ) ; % number of rows of the matrix
Dxy = image_width / ( n - 1 ) ;
dt = 1e-12 ; % differential of time to perform derivation
tpdt = t + dt ; % time instants for t+dt
t = t - dt ;

sparsityFactor = .011 * n / 150 ;
A_mat = sparse( [], [], [], nn, n_rows, round( nn * n_rows * sparsityFactor ) ) ;
% A_mat = [] ;

dist_sensor = c0 * t ;
dist_pdd_sensor = c0 * tpdt ;

% progress bar
wbar = waitbar(0,'','Name','Calculating model matrix','CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');

r_counter = 1 ;
perproj = 0;
for ii = 1 : length( angle_sensor )
    
    if getappdata(wbar,'canceling')
        A_mat = [];
        Amatfile = [];
        close(wbar);
        delete(wbar);
        error('User cancelled the creation of model matrix');
        return;
    end
    
    % time estimation
    est = round(perproj*(length( angle_sensor )-ii+1)); est_unit = 's';
    if (est > 120) est = round(est / 60); est_unit = 'min'; end;
    
    
    wbar = waitbar(ii/length( angle_sensor ),wbar,sprintf( 'Estimated Time Left: %i%s', est, est_unit )) ;
    tic;
    
    angle_max = asin( ( ( image_width + 2 * Dxy ) * sqrt(2) ) / ( 2 * r_sensor( r_counter ) ) ) ;
    if( ~isreal( angle_max ) )
        angle_max = pi ;
    end ;
    angles = linspace( -angle_max, angle_max, n_angles )' * ones( 1, length( t ) ) ;
        
    theta = -angle_sensor( ii ) ;
    
    x_pt = r_sensor( r_counter ) - ( ones( n_angles, 1 ) * dist_sensor ) .* cos( angles ) ;
    y_pt = ( ones( n_angles, 1 ) * dist_sensor ) .* sin( angles ) ;
    R_pt = ones( n_angles, 1 ) * dist_sensor ;
    xpdx_pt = r_sensor( r_counter ) - ( ones( n_angles, 1 ) * dist_pdd_sensor ) .* cos( angles ) ;
    ypdy_pt = ( ones( n_angles, 1 ) * dist_pdd_sensor ) .* sin( angles ) ;
    RpdR_pt = ones( n_angles, 1 ) * dist_pdd_sensor ;
    
    A_tmp = ( 1 / ( 2*dt ) ) * ( -computeDiscretizedMatrix_proj( x_pt, y_pt, R_pt, theta, image_width, n, lt ) + computeDiscretizedMatrix_proj( xpdx_pt, ypdy_pt, RpdR_pt, theta, image_width, n, lt ) ) ;
    A_mat( :, (ii-1)*lt+1:ii*lt ) = A_tmp ;
    
    if( length( r_sensor ) > 1 )
        r_counter = r_counter + 1 ;
    end ;
    perproj = toc;
end


close(wbar);
delete(wbar);










