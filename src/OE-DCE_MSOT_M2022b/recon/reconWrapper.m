function [ Recon, A_mat, sigMat_filtered, return_code, return_message ] = reconWrapper( A_mat, sigMat, imp_resp, n, proj, r_sensor, angle_sensor, c, image_select, fs, limits, time_res, filter_f, lsqr_iter, varargin )
%% INPUT
% @param A_mat: forward model-based matrix (sparse and possibly convolved,
%               then imp_resp should be 0).
%               Or inverted (recon) matrix in full (not sparse).
% @param sigMat: matrix of measured signals, each transducer one column.
%                Time starts at row 1 and ends at last row.
% @param imp_resp: impulse response of the transducer, 0 or [] if N/A.
% @param n: pixels of square image in one dimension
% @param proj: projections
% @param r_sensor: the distance of sensor(s) from the center of the coord.
%                  Can be a vector, if every sensor has a different
%                  distance. in [m]
%                  Can be Matrix, if cartesian coordinates are used. 
%                  In this case size is [proj x 3].
% @param angle_sensor: the angle of the sensor(s) relative to the
%                      coordinate system. in [RAD]
%                      Can be empty, if cartesian coordinates are used.
% @param c: average speed of sound in sample. in [m/s]
% @param image_select: select the reconstructed image. Values are 'full', 
%                      'direct', 'derivative' or 'model_lin'.
% @param fs: sampling frequency. in [Hz]
% @param limits: a vector with two entries [start end], describing the start/end of the acoustic signal. in [m]
% @param time_res: time resolution for model-based reconstruction.
% @param filter_f: two element vector containing cutt-off frequencies for
%                  filter_f(1) - low-pass filter
%                  filter_f(2) - high-pass filter
%                  optional: supply third entry (1) for non-zero-phase
%                            filtering
% @param lsqr_iter: number of LSQR iterations. Choose very low value for
%                   limited view geometry to prevent ripples.
% @opt_param pwrCorr: power correction factor. sigMat will be devided by this.
% @opt_param coreNum: number of cores to use
% @opt_param verbose: turn on/off diagnostic messages.
% @opt_param progressBar: JProgressBar to show progress in a GUI.
%
% @return Recon: reconstructed image n-by-n.
% @return A_mat: same as input A_mat or new if empty array was passed. If an
%                impulse response was also passed, A_mat is also convolved.
% @return sigMat_filtered: uncut, filtered signal matrix.
% @return return_code, equals 0 if no error. -1 if warning and -2 if error.
% @return_message, message of the exception if any.
%
% Temp to Speed: c = 1.402385 * 1e3 + 5.038813 * T - 5.799136 * 1e-2 * T^2 + 3.287156 * 1e-4 * T^3 - 1.398845 * 1e-6 * T^4 + 2.787860 * 1e-9 * T^5 ;

% persistent sigMat_figure ;

persistent filter_struct ;

if( isempty( filter_struct ) )
    
    filter_struct = struct( 'passband', [0 0], ...
                                    'a_HPF', [], ...
                                    'b_HPF', [], ...
                                    'a_LPF', [], ...
                                    'b_LPF', [] ) ;
end ;

try
    %%%%%%%%%%%%%%%%%%%%%%%% Return params init %%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigMat_filtered = [] ;
    Recon = eye(n) ;
    return_code = 0 ;
    return_message = '' ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pwrCorr = 1 ;
    coreNum = feature( 'numCores' ) ;
    verbose = 0 ;
    progressBar = [] ;
    
    if( nargin > 14 )
        vlen = length( varargin ) ;
        if( vlen > 0 )
            pwrCorr = varargin{1} ;
            if( pwrCorr == 0 )
                pwrCorr = 1 ;
            end ;
        end ;
        if( vlen > 1 )
            coreNum = varargin{2} ;
        end ;
        if( vlen > 2 )
            verbose = varargin{3} ;
        end ;
        if( vlen > 3 )
            progressBar = varargin{4} ;
        end ;
    end ;
    
    if( verbose )
        tic
    end ;
    
    sigMat = sigMat / pwrCorr ;

    %% ------------------------ PARAMETER CHECK      -------------------------%
    if( ~sum( strcmp(image_select, {'full' 'direct' 'derivative' 'model_lin' 'wavelet'} ) ) )
        warning('RECONWRAPPER:reconValueNotAllowed', '[RECONWRAPPER]; No valid reconstruction mode selected, using full-backprojection.') ;
        image_select = 'full' ;
        time_res = 10 ;
        return_message = [ return_message 'No valid reconstruction mode selected, using full-backprojection;' ] ;
    end ;

%     if( isvector( A_mat ) )
%         error( 'RECONWRAPPER:wrongReconMatrix', 'Reconstruction matrix must not be a number or a vector.' ) ;
%     end ;

    if( size( sigMat, 1 ) < size( sigMat, 2 ) )
       sigMat = sigMat' ;
    end ;

    if( length( angle_sensor ) ~= size( sigMat, 2 ) && size( r_sensor, 1 ) ~= size( sigMat, 2 ) )
        error( 'RECONWRAPPER:dataMismatch', 'Mismatch between sigMat (data) and description of sensor positions.' ) ;
    end ;

    if( time_res <= 1 )
        warning( 'RECONWRAPPER:valueDefaulted', '[RECONWRAPPER]; Time resolution found to be smaller than 1. Reset to 3.' ) ;
        time_res = 3 ;
        return_message = [ return_message 'Time resolution found to be smaller than 1. Reset to 3;' ] ;
    end ;

    if( max( size( imp_resp ) ) > 1 )
        if( size( imp_resp, 1 ) < size( imp_resp, 2 ) )
            imp_resp = imp_resp' ;
        end ;
    end ;

    % ------------------------- PARAMETER CHECK END --------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ------------------------ DECONVOLUTION --------------------------------%
    if ( max( size( imp_resp ) ) > 1 )
        if( verbose )
            basetime = toc ;
        end ;
        
        firstNNZ = find( sigMat( :, 1 ) ~= 0, 1 ) ; firstNNZ = firstNNZ + 10 ;
        sigMat( 1:firstNNZ-1, : ) = 0 ;
        sigMat( firstNNZ:end, : ) = sigMat( firstNNZ:end, : ) - ones( length( sigMat( firstNNZ:end, 1 ) ), 1 ) * sigMat( firstNNZ, : ) ;
        
        for kk = 1 : size( sigMat, 2 )
            sigMat( :, kk) = deconvwnr( sigMat( :, kk ), imp_resp, 0 ) ;
        end ;
        
        if( verbose )
            time = toc - basetime ;
            fprintf( ['\n[RECONWRAPPER]; Deconvolution took ' num2str( time*1e3 ) ' ms. \n'] ) ;
        end ;
        
    end ;
    % ------------------------- DECONVOLUTION END ----------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ------------------------ Pre-Filter Optimization ----------------------%
    % remove trigger/backscattering artefact from the first samples
    if( verbose )
    	basetime = toc ;
    end ;
        
%     firstNNZ = find( sigMat( :, 1 ) ~= 0, 1 ) ; 
%     if( ~isempty( firstNNZ ) && firstNNZ < 50 )
    firstNNZ = 150 ;
    sigMat( 1:firstNNZ-1, : ) = 0 ;
    sigMat( firstNNZ:end, : ) = sigMat( firstNNZ:end, : ) - ones( length( sigMat( firstNNZ:end, 1 ) ), 1 ) * sigMat( firstNNZ, : ) ;
%     end ;
    if( filter_f(2) ~= 0 || filter_f(1) ~= 0 )
        fsi = 3 * fs ;
        t = 1/fs : 1/fs : size( sigMat, 1 ) / fs ;
        ti = t(1) : 1/fsi : t(end) ;
        t_norm = ( ti - t(1) ) * fs + 1 ;
        [XI YI] = meshgrid( 1 : size( sigMat, 2 ), t_norm ) ;
        sigMat = ba_interp2( sigMat, XI, YI, 3, coreNum ) ;
        fs = fsi ;
        clear XI YI X Y t ti t_norm ;
    end ;
    
    if( verbose )
        time = toc - basetime ;
        fprintf( ['\n[RECONWRAPPER]; Upsampling took ' num2str( time*1e3 ) ' ms. \n'] ) ;
    end ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% --------------------------- SOME SETUP --------------------------------%
    image_width = limits(2) - limits(1) ;
    dx = image_width/n ;
    dt = dx/(time_res*c) ;

    if( filter_f(2) ~= 0 )
        if( filter_struct.passband(2) ~= filter_f(2) || isempty( filter_struct.a_LPF ) )
            f_LPF = filter_f(2) ;
            [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .9 ) ;
            filter_struct.passband(2) = filter_f(2) ;
            filter_struct.a_LPF = a_LPF ;
            filter_struct.b_LPF = b_LPF ;
        else
            b_LPF = filter_struct.b_LPF ;
            a_LPF = filter_struct.a_LPF ;
        end ;
    end ;

    if( filter_f(1) ~= 0 )
        if( filter_struct.passband(1) ~= filter_f(1) || isempty( filter_struct.a_HPF ) )
            f_HPF = filter_f(1) ;
            [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;
            filter_struct.passband(1) = filter_f(1) ;
            filter_struct.a_HPF = a_HPF ;
            filter_struct.b_HPF = b_HPF ;
        else
            b_HPF = filter_struct.b_HPF ;
            a_HPF = filter_struct.a_HPF ;
        end ;
    end ;

    len = size(sigMat, 1) ;
    ts = 1/fs:1/fs:len/fs ;

    fac = fs/c ;
    sig_start = int32( (limits(1) ) * fac - (sqrt(2)-1)*image_width/2 * fac ) ;
    sig_end = int32( (limits(2) ) * fac + (sqrt(2)-1)*image_width/2 * fac ) ;

    if( sig_start == 0 )
        sig_start = 1 ;
    end;
    % ---------------------------- END SETUP ---------------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ------------------------- FILTER --------------------------------------%
    
    if( verbose )
        basetime = toc ;
    end ;
    
    if( filter_f(2) ~= 0 )
        if (numel(filter_f) > 2)
            sigMat = FilterM(b_LPF,a_LPF,sigMat);
        else
            sigMat = FiltFiltM( b_LPF, a_LPF, sigMat, 1, coreNum ) ;
        end
        cut = ceil( .05 * size( sigMat, 1 ) ) ;
        sigMat = sigMat - ones( size( sigMat, 1 ), 1 ) * sigMat( cut+1, : ) ;
        sigMat( [1:cut end-cut:end], : ) = 0 ;
    end;

    if( filter_f(1) ~= 0 )
        if (numel(filter_f) > 2)
            sigMat = FilterM(b_HPF,a_HPF,sigMat);
        else
            sigMat = FiltFiltM( b_HPF, a_HPF, sigMat, 1, coreNum ) ;
        end
    end;
    
    sigMat_filtered = sigMat ;
    
    if( verbose )
        time = toc - basetime ;
        fprintf( ['\n[RECONWRAPPER]; Filtering took ' num2str( time*1e3 ) ' ms. \n'] ) ;
    end ;

    % -------------------------- FILTER/MEAN/NORMALIZE END -------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ------------------------- CUT AND DOWNSAMPLE --------------------------%
    % cut signal and impulse response (if available)
    
    if( verbose )
        basetime = toc ;
    end ;
    
    % set up AA-Lowpass
    r = ( ts(2)-ts(1) ) / dt * 3 ; % corresponds to ( f_ds / fs )
    if( r <= 1 )
        [b,a] = cheby1(8, .05, .938*r) ;
        upsample = 0 ;
    else
        upsample = 1 ;
    end ;

    % low-pass filter
    if( ~upsample )
        sigMat = FiltFiltM( b, a, sigMat, 1, coreNum ) ;
    end ;
    
    num_proj = size( sigMat, 2 ) ;
    if( max( size( r_sensor ) ) > 1 && max( size( r_sensor ) ) == num_proj && ~isempty( angle_sensor ) )

        % LINESCAN : first resample
        if ( num_proj ~= proj )

            from = r_sensor(1) * cos( angle_sensor(1) ) ;
            to = r_sensor(end) * cos( angle_sensor(end) ) ;
            data_pos = linspace( from, to, num_proj ) ;
            intp_pos = linspace( from, to, proj ) ;

            r_intp = interp1( data_pos, r_sensor, intp_pos, 'cubic' ) ;
            phi_intp = interp1( data_pos, angle_sensor, intp_pos, 'cubic' ) ;

        else
            phi_intp = angle_sensor ;
            r_intp = r_sensor ;
        end ;

        t = ts(1) : dt : ts(end) ;

        [XI YI] = meshgrid( phi_intp, t ) ;
        [X Y] = meshgrid( angle_sensor, ts ) ;
        sigMat2 = interp2( X, Y, sigMat, XI, YI, 'cubic' ) ;
        angle_sensor = phi_intp ;
        r_sensor = r_intp ;

        % LINESCAN : then cut
        sig_start = round( limits(1) / c / dt ) ;
        sig_end = round( limits(2) / c / dt ) ;

        if( sig_start < 1 )
            sig_start = 1 ;
        end ;
        if( sig_end > size( sigMat2, 1 ) )
            sig_end = size( sigMat2, 1 ) ;
        end ;

        sigMat2 = sigMat2( sig_start:sig_end, : ) ;
        t = t( sig_start:sig_end ) ;

    else
        try

            if( sig_end > size( sigMat, 1 ) )
                sig_end = size( sigMat, 1 ) ;
            end ;

            if( sig_start < 1 )
                sig_start = 1 ;
            end ;

            sigMat = sigMat( sig_start:sig_end, : ) ;
            ts = ts( sig_start : sig_end ) ;

        catch ex
            ex = addCause( ex, MException( 'RECONWRAPPER:ReconLimitsException', '[RECONWRAPPER]; Recon limits too big. Please adjust region of interest.' ) ) ;
            rethrow( ex ) ;
        end ;
    end ;

    % re-sample and interpolate if necessary
    if( ~isempty( angle_sensor ) && max( size( r_sensor ) ) == 1 )
        
        % polar coordinates, i.e. Tomography
        t = ts(1) : dt : ts(end) ;

        phi_intp = angle_sensor ;
        if ( num_proj ~= proj )
            dphi =  (angle_sensor(end)-angle_sensor(1) ) / (proj-1) ;
            phi_intp = angle_sensor(1) : dphi : angle_sensor(end) ;
        end ;

        if( length( angle_sensor ) > 1 )
            phi_intp_norm = ( phi_intp - angle_sensor(1) ) / (angle_sensor(2) - angle_sensor(1)) + 1 ;
        else
            phi_intp_norm = 1 ;
        end ;

        t_norm = ( t - ts(1) ) * fs + 1 ;
        [X Y] = meshgrid( phi_intp_norm, t_norm ) ;
        sigMat2 = ba_interp2( sigMat, X, Y, 3, coreNum ) ;
%         [XI YI] = meshgrid( phi_intp, t ) ;
%         [X Y] = meshgrid( angle_sensor, ts ) ;        
%         sigMat2 = interp2( X, Y, sigMat, XI, YI, 'spline' ) ;
        angle_sensor = phi_intp ;

    elseif( isempty( angle_sensor ) && size( r_sensor, 2 ) == 3 )
        
        % cartesian coordinates in 3D, i.e. free geometry (no projections re-sampling)
        t = ts(1) : dt : ts(end) ;
        phi_intp_norm = 1 : size( r_sensor, 1 ) ;
        t_norm = ( t - ts(1) ) * fs + 1 ;
        [X Y] = meshgrid( phi_intp_norm, t_norm ) ;
        sigMat2 = ba_interp2( sigMat, X, Y, 3, coreNum ) ;

    end ;

    % local clean-up
    clear dphi phi_intp t_norm X Y ;
    clear X Y n_intp start_index a b r dx ts ;
    clear a_LPF b_LPF a_HPF b_HPF f_LPF f_HPF fac ii len num_proj ;
    clear sig_end sig_start ;
    
    if( verbose )
        time = toc - basetime ;
        fprintf( ['\n[RECONWRAPPER]; Cut&Downsample took ' num2str( time*1e3 ) ' ms. \n'] ) ;
    end ;
    % ------------------------ CUT AND DOWNSAMPLE END ------------------------%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ------------------------------ DO RECON -------------------------------%
    if ( strcmp( 'model_lin', image_select ) )

        % do model-based recon
        if( isempty(A_mat) )
            
            if( ~isempty( progressBar ) )
                progressBar.setString( 'Building forward matrix.' ) ;
            end ;
            
%             A_mat = calcTriAngles( c, n, image_width, t, r_sensor, angle_sensor, progressBar ) ;
            
            A_mat = computeDiscretizedMatrix( c, n, image_width, t, r_sensor, angle_sensor, 2*n ) ;
            
        else
            
            [rows cols] = size( sigMat2 ) ;
            equations = max( size( A_mat ) ) ;

            if( cols*rows ~= equations )
                offset =  equations/cols - rows ;
                if( (rows + round(offset) )*cols ~= equations || abs( offset ) > 5 )
                    error( [ 'RECONWRAPPER:matrixError', 'Matrix does not fit data. Time sample offset is ' num2str( offset ) '. Reconstruction aborted ...' ] ) ;
                else
                    if( offset < 0 )
                        sigMat2( rows+offset+1:rows, : ) = [] ;
                    elseif( offset > 0 )
                        sigMat2( rows+offset, : ) = zeros( offset, cols ) ;
                    end ;
                end ;
            end ;
            
            [rows cols] = size( sigMat2 ) ;
            if( cols*rows ~= equations )
                error( [ 'RECONWRAPPER:matrixError', 'Matrix does not fit data. Correction did not work: cols = ' num2str( cols ) '; rows = ' num2str( rows ) '; equations = ' num2str( equations ) ';' ] ) ;
            end ;

        end ;
        
        [userview systemview] = memory ;
        maxMem = systemview.PhysicalMemory.Total ;
        whosAmat = whos( 'A_mat' ) ;
        if( whosAmat.bytes < maxMem && issparse( A_mat ) )
            %% SPARSE matrices are expected to be FORWARD matrices
            if( size( A_mat, 1 ) > size( A_mat, 2 ) )
                
                if( verbose )
                    fprintf('\n[RECONWRAPPER]; Starting LSQR with max %d iterations ... \n\n', lsqr_iter) ;
                end ;
                
                [msg_lsqr recon] = evalc( 'lsqr( A_mat, ( reshape(sigMat2, size(sigMat2, 1)*size(sigMat2, 2), 1) ), 1e-6, lsqr_iter)' ) ;
                
            else
                
                if( verbose )
                    fprintf( '\n[RECONWRAPPER]; Starting LSQR_mod with max %d iterations ... \n\n', lsqr_iter ) ;
                end ;
                
                setappdata( 0, 'Recon_mat', A_mat ) ;
                setappdata( 0, 'LSQR_progress', struct( 'isGUI', ~isempty( progressBar ), 'iter', -1 ) ) ;
                [msg_lsqr recon] = evalc( 'lsqr( @lsqr_mod, ( reshape(sigMat2, size(sigMat2, 1)*size(sigMat2, 2), 1) ), 1e-6, lsqr_iter)' ) ;
                rmappdata( 0, 'Recon_mat' ) ;
                rmappdata( 0, 'LSQR_progress' ) ;
                
            end ;

            if( verbose )
                fprintf( ['\n[RECONWRAPPER]; ' msg_lsqr '\n'] ) ;
            end ;
            
            Recon = reshape(recon, n, n) ;

        end ;

    % ****************************************************************
    % Wavelet Reconstruction, SM, 23/01/2013
    % ****************************************************************
    elseif ( strcmp( 'wavelet', image_select ) )
        % if no input, first set up parameters
        if (isempty(A_mat))
            A_mat = struct();
            A_mat.thres = 0.2;      % threshold to ignore crosstalk betw WL
            A_mat.over_fact = 1.4;  % select timepoints for each WL
            A_mat.Athres = 0.01;    % thresholding of inv mat to be sparse
            A_mat.Sthres = 0;       % per default, do not threshold SVD
            A_mat.depth = 2;        % default is two levels of decomp.
            A_mat.wl_name = 'db6';  % name of the mother wavelet
        end
        
        % make signal vector 1D
        sigMat2 = reshape(sigMat2,size(sigMat2,1)*size(sigMat2,2),1);
        
        % convert signal to wl domain
        p_w = transform_t2wl( sigMat2, numel(t),numel(angle_sensor),A_mat.depth,A_mat.wl_name);
      
        
        % if no inverted set of matrices given, calculate first the model
        % and then the inverse in the WL domain
        if (~isfield(A_mat,'Ri') || isempty(A_mat.Ri))
            % calculate forward model
            FM = computeDiscretizedMatrix( c, n, image_width, t, r_sensor, angle_sensor, 2*n ) ;
            FM = FM';
            A_mat = WL_invert(A_mat,p_w,FM,n,t,numel(angle_sensor),verbose);
            clear FM;       % delete forward model
        end
        
        % reconstruct image
        if (verbose)
            fprintf('Starting wl-domain reconstruction...\n');
        end
        Recon=zeros(n,n);
        
        for kk=1:4^(A_mat.depth)    
            if (verbose) fprintf('Wavelet packet %i of %i\n',kk,4^(A_mat.depth)); end;
        
            % use relevant places only
            pi_w_dash=p_w(A_mat.places{kk}).';

            % reconstruction in wavelet domain, multiply inverted model matrix with
            % measured signal
            recon_wlet=A_mat.Mi_w_dash_inv{kk}*sparse(pi_w_dash(:));

            % transform to spatial domain and add to previous wavelets
            Recon = Recon + reshape( A_mat.Ri{kk}*recon_wlet(:),n,n);
        end
        clear kk recon_wlet;
        
    else
        %% do backprojection
        image_select_bp = 3 ;
        if( strcmp( 'direct', image_select ) )
            image_select_bp = 1 ;
        elseif( strcmp( 'derivative', image_select ) )
            image_select_bp = 2 ;
        end;

        if( size( r_sensor, 2 ) ~= 3 )
            if( size( angle_sensor, 2 ) ~= 1 )
                r_sensor = r_sensor' ;
                angle_sensor = angle_sensor' ;
            end ;
            if( length( r_sensor ) == 1 )
                r_sensor = [ cos( angle_sensor ) sin( angle_sensor ) zeros( size( angle_sensor ) ) ] * r_sensor ;
            else
                r_sensor = [ cos( angle_sensor ).* r_sensor sin( angle_sensor ).* r_sensor zeros( size( angle_sensor ) ) ] ;
            end ;
        end ;

        if( size( sigMat2, 1 ) ~= size( t, 1 ) )
            t = t' ;
        end ;
        
        if( verbose )
            basetime = toc ;
        end ;
        
        Recon = - backproject( sigMat2, n, r_sensor(:, 1), r_sensor(:, 2), r_sensor(:, 3), c, image_select_bp, t, 1/dt, image_width, coreNum ) ;
        
        if( verbose )
            time = toc - basetime ;
            fprintf( ['\n[RECONWRAPPER]; BP took ' num2str( time*1e3 ) ' ms. \n'] ) ;
        end ;
    end ;
    
    return_message = [return_message 'Reconstruction successful.' ] ;
    
catch ex
    
    clear sigMat ;
    return_code = -2 ;
        
    msg = ex.getReport ;
    msg = regexprep( msg, '</a>', '' ) ;
    msg = regexprep( msg, '<a href.*>', '' ) ;
    
    return_message = [ '[RECONWRAPPER-ERROR-REPORT] ' char( msg ) ] ;
        
    if( ~isdeployed )
        rethrow(ex) ;
    end ;
    
end ;







