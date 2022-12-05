function y = lsqr_mod( x, transp_flag )

   A = getappdata( 0, 'Recon_mat' ) ;
   
   if strcmp(transp_flag,'transp')      % y = A'*x 
       
       lsqr_progress = getappdata( 0, 'LSQR_progress' ) ;
   
       lsqr_progress.iter = lsqr_progress.iter + 1 ;
       if( lsqr_progress.iter )
%            if( lsqr_progress.isGUI )
%                java.lang.System.out.patprintln( java.lang.String( [ 'Starting iteration ' num2str( lsqr_progress.iter ) ] ) ) ;
%                drawnow ;
%            else
               fprintf( [ '\nStarting iteration ' num2str( lsqr_progress.iter ) '\n' ] ) ;
               drawnow ;
%            end ;
       end ;
       setappdata( 0, 'LSQR_progress', lsqr_progress ) ;
       
       y = A * x ;
       
   elseif strcmp(transp_flag,'notransp') % y = A*x -> (A*x)' = x'*A'
       
       if( size(x, 2) == 1 )
           y = x' * A ;
       else
           y = x * A ;
       end ;
   end
   
   if( size( y, 2 ) ~= 1 )
       y = y' ;
   end ;
   
end