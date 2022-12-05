function [ image_filt ] = filterImage3D( varargin )
%% Filters image frequencies by 3D hanning windowing of the frequency transform.
%
%  @param image, the image to filter
%  @param filter_f, the 2 filter coefficients for high/low-pass filtering.
%                   [ High_F Low_F ] => Bandpass
%
%                   For HIGH-PASS: the higher the coefficient, the lower
%                   the cut-off frequency.
%                   For LOW-PASS: the lower the coefficient, the higher the
%                   cut-off frequency.
%
%  @return image_filt, the filtered image.
%
%  usage: image_filt = filterImage3D( image, [1 2] ) ;
%  self-test: filterImage() ;

if( nargin == 0 )
    selftest() ;
    image_filt = 0 ;
    return
end ;

image = varargin{1} ;
filter_f = varargin{2} ;
offset= varargin{3};

[n m k] = size( image ) ;
P_f0 = fftshift( fftn( ifftshift( image + 0 ) ) );

hanning3d = getWin([n m k], 'Hanning') ;

if( filter_f(2) ~= 0 && filter_f(1) == 0 )
    
    P_f0 = P_f0 .*hanning3d.^filter_f(2) ;

elseif( filter_f(1) ~= 0 && filter_f(2) == 0 )
    
    P_f0 = P_f0 .* ( 1 - hanning3d.^filter_f(1) ) ;

elseif( filter_f(1) ~= 0 && filter_f(2) ~= 0 )
    
    LP = hanning3d.^filter_f(2) ;
    HP = 1 - hanning3d.^filter_f(1) ;
    weight = LP .* HP ;
    P_f0 = P_f0 .* weight ;

else
    image_filt = image ;
    return 
end ;

image_filt = real( fftshift( ifftn( ifftshift( P_f0 ) ) ) ) ;
T=zeros(size(image_filt));
T(1+offset:end-offset,1+offset:end-offset,1+offset:end-offset)=1;
image_filt=T.*image_filt;

end

function selftest()
    k = 244;
    r_ph = 0.0125 ;
    dx = 0.03/(200-1);
    x = -0.03/2:dx:0.03/2;
    y = x;
    [X Y] = meshgrid(x,y);
    r = sqrt(X.^2+Y.^2);
    U = besseli(0, k*r).*( r < r_ph );
    
    figure( 'Name', 'Original' ) ; imagesc( U ) ; colorbar ;
    mainf = figure( 'Name', 'Demo Usage', 'Units', 'Normalized' ) ;
    pos = get( mainf, 'Position' ) ;
    pos2 = get( mainf, 'OuterPosition') ;
    set( mainf, 'Position', [ 0 pos(2)-pos2(4) 1 pos(4) ] ) ;
    counter = 1 ;
    for ii = [ .1 .5 1 5 10 50 100 ]
        
        h = subplot( 2, 7, counter ) ;
        imagesc( filterImage( U, [ ii 0 ] ) ) ; colorbar ; axis image ;
        title( [ 'HP ' num2str(ii) ], 'FontSize', 10, 'FontWeight', 'bold' ) ;
        
        h = subplot( 2, 7, 7 + counter ) ;
        imagesc( filterImage( U, [ 0 ii ] ) ) ; colorbar ; axis image ;
        title( [ 'LP ' num2str(ii) ], 'FontSize', 10, 'FontWeight', 'bold' ) ;
        
        counter = counter + 1 ;
    end ;
end














