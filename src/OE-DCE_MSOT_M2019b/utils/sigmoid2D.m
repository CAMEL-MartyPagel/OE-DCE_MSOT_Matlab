function imnew = sigmoid2D(im, a, b)
% Apply sigmoid intensity transformation

% 
%       ARGUMENT DESCRIPTION:
%               IM       - gray scale image (MxN).
%               a, b     - parameters
% 
%       OUTPUT DESCRIPTION:
%                IMNEW - output image .
% 
   Imax = max(im(:));
   Imin = min(im(:));
   imnew = (Imax - Imin)*(1./(1+exp(-((im-b)/a)))) + Imin;
   imnew = imnew.^(2.2);
   imnew = newrange(imnew,[Imin Imax]);

end