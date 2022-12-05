function sigMatf = filterData( sigMat_full, filter_f, fs, varargin )

interp = true;
if (numel(varargin)>=1 && ~varargin{1})
    interp = false;
end

f_LPF = filter_f(2) ;
f_HPF = filter_f(1) ;
if interp,
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs/3 * .9 ) ;
    [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs/3 * 1.46, 'high' ) ;
else
    [b_LPF,a_LPF] = cheby1( 8, .01, 2 * f_LPF/fs * .9 ) ;
    [b_HPF,a_HPF] = cheby1( 4, .01, 2 * f_HPF/fs * 1.46, 'high' ) ;
end    

% fprintf( '\n' ) ;
svec = size(sigMat_full);
sigMat_full = reshape(sigMat_full,[svec(1) svec(2) prod(svec(3:end))]);

if interp,
    sigMatf = zeros(6088,size(sigMat_full,2),size(sigMat_full,3));
else
    sigMatf = zeros(svec(1),svec(2),size(sigMat_full,3));
end


for ii = 1 : size( sigMat_full, 3 )
    
%     fprintf( [ 'Filtering ' num2str(ii) ' of ' num2str( size( sigMat_full, 3 ) ) '\n' ] ) ;
    
    sigMat = sigMat_full( :, :, ii ) ;
    
    if interp
        firstNNZ = find( sigMat( :, 1 ) ~= 0, 1 ) ; firstNNZ = firstNNZ + 10 ;
        sigMat( 1:firstNNZ-1, : ) = 0 ;
        sigMat( firstNNZ:end, : ) = sigMat( firstNNZ:end, : ) - ones( length( sigMat( firstNNZ:end, 1 ) ), 1 ) * sigMat( firstNNZ, : ) ;
        fsi = 3 * fs ;
        t = 1/fs : 1/fs : size( sigMat, 1 ) / fs ;
        ti = t(1) : 1/fsi : t(end) ;
        t_norm = ( ti - t(1) ) * fs + 1 ;
        [XI YI] = meshgrid( 1 : size( sigMat, 2 ), t_norm ) ;
        sigMat = ba_interp2( sigMat, XI, YI, 3, 2 ) ;
        fs = fsi ;
        clear XI YI X Y t ti t_norm ;
    end
    
    sigMat = FiltFiltM( b_LPF, a_LPF, sigMat, 1, 2 ) ;
    
    if interp
        cut = ceil( .05 * size( sigMat, 1 ) ) ;
        sigMat = sigMat - ones( size( sigMat, 1 ), 1 ) * sigMat( cut+1, : ) ;
        sigMat( [1:cut end-cut:end], : ) = 0 ;
    end
    
    sigMatf( :, :, ii ) = FiltFiltM( b_HPF, a_HPF, sigMat, 1, 2 ) ;
    
end ;
    
if interp,
    sigMatf = reshape(sigMatf,[6088 svec(2:end)]);
else
    sigMatf = reshape(sigMatf,svec);
end
