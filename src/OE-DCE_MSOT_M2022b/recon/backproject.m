function Recon = backproject( sigMat, n, x, y, z, c, image_select_bp, t, fs, image_width, coreNum )
%% MEX implementation of Backprojection, i.e. parallelized machine code.
%  All variables in SI units.
%  
% @param sigMat, data in samples-by-channels
% @param n, pixels of square image along one dimension
% @param x, x-values of each channel in cart. coord. Size channels-by-1
% @param y, y-values of each channel in cart. coord. Size channels-by-1
% @param z, z-values of each channel in cart. coord. Size channels-by-1
% @param c, speed of sound
% @param, image_select_bp, 1 direct, 2 derivative, 3 full
% @param, t, time vector describing the first dimension of sigMat
% @param, fs, sampling rate of sigMat
% @param, image_width, width of the reconstructed image around the origin.
% @param, coreNum, number of cores to use. If set to more than available the code will slow down.
% 
% @return, Recon, reconstructed image.


% function [ Recon ] = backproject( sigMat, n, r_sensor, angle_sensor, c, image_select, t, fs, image_width, numCores )

% if( isempty( strmatch(image_select, {'full' 'direct' 'derivative'}, 'exact') ) )
%     warning('MATLAB:valueNotAllowed', 'No valid backprojection image mode selected, using full.');
%     image_select = 'full';
% end;

% sizeT = length(t);
% 
% bp_image = zeros(n, n);
% bp_image2 = zeros(n, n);
% bp_image_12 = zeros(n, n);
% 
% [N, M] = meshgrid( 1:n, 1:n );
% xr = ( N - (n+1)/2 ) * ( image_width/(n-1) ) ;
% yr = ( (n+1)/2 - M ) * ( image_width/(n-1) ) ;
% 
% for i = 1:length(angle_sensor)
%     
%     A0 = sigMat(:,i);
%     phi = angle_sensor(i);
%     
%     if ( length(r_sensor) == 1 )
%         x_sensor = r_sensor ;
%     else
%         x_sensor = r_sensor(i) ;
%     end;
%     y_sensor = 0 ;
%     
%     x = x_sensor - (xr*cos(-phi) - yr*sin(-phi));
%     y = y_sensor - (xr*sin(-phi) + yr*cos(-phi));
%     s = sqrt( x.^2 + y.^2 );
%     ss = round( s*fs/c - t(1)*fs + 1 );
%     ss( ss > sizeT ) = sizeT;
%     ss( ss <= 0 ) = 1;
%     
%     s1 = A0(ss);
% %     A4 = [A0(2:sizeT)'-A0(1:sizeT-1)' 0]'./t';
%     A4 = [A0(2:sizeT)'-A0(1:sizeT-1)' 0]'.*t'/(t(2)-t(1));
%     s2 = A4(ss);
%     bp_image = bp_image + s1;
%     bp_image2 = bp_image2 + s2;
%     bp_image_12 = bp_image_12 + s1 - s2;
% end ;
% 
% if( image_select == 1 )
%     Recon = bp_image;
% elseif( image_select == 2 )
%     Recon = bp_image2;
% else
%     Recon = bp_image_12;
% end;
















