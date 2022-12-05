function [x, f] = mbd(x, y, maxiter, varargin);
%MBD does a single for Multi-Frame Blind Deconvolution.
%
% Purpose:
%   Implements the algorithm in:
%
%      Stefan Harmeling, Michael Hirsch, Suvrit Sra, Bernhard Sch\"olkopf,
%      "Online Blind Image Deconvolution for Astronomy", in Proceedings of
%      the First International Conference on Computational Photography, 2009.
%
% Inputs:
%   x              - current estimate of the deblurred image
%                    possible initialization is ones(sx)
%   y              - next observed image
%   maxiter        - e.g. [100, 3] means 100 update steps for f and 3 for x
%
% Note that the size of x determines the size of the PSF, via
%         size(f) = size(x) - size(y) + 1
%
% Outputs:
%   x              - new estimate of the deblurred image
%
% Copyright (C) 2009 by Stefan Harmeling (2009-09-23).

sx = size(x);            % size of PSF             
sy = size(y);            % size of blurred image   
sf = sx - sy + 1;        % size of true image      

if numel(varargin)>0 && ~isempty(varargin{1})
   f = varargin{1};
else
% estimate PSF f
f = norm(y(:)) / norm(x(:));
f = f * ones(sf) / sqrt(prod(sf));
tic
f = mbd_update(f, x, y, maxiter(1));   % estimate the PSF
toc
sumf = sum(f(:));
f = f/sumf;                         % normalize f
x = sumf*x;                         % adjust x as well
end

% update true image x
x = mbd_update(x, f, y, maxiter(2));   % estimate the image
return

%%%%%%%%%%%%%%%%
function f = mbd_update(f, x, y, maxiter)
% depending on the value of sf the roles of f and x can be swapped
sy = size(y);
sf = size(f);
for i = 1:maxiter
  ytmp = pos(cnv2(x, f));
  nom = pos(cnv2tp(x, y));
  denom = pos(cnv2tp(x, ytmp));
  tol = 1e-10;
  factor = nom ./(denom+tol);
  sel = find(abs(nom-denom) < tol);
  factor(sel) = 1;
  factor = reshape(factor, sf);
  f = f .* factor;
end
return

%%%%%%%%%%%%%%%%%
function x = pos(x, epsilon)
x(find(x(:)<0)) = 0;
return

%%%%%%%%%%%%%%%%%
function A = cnv2slice(A, i, j);
A = A(i,j);
return

%%%%%%%%%%%%%%%%%
function y = cnv2(x, f)
sx = size(x);
sf = size(f);
if all(sx >= sf)   % x is larger or equal to f
  % perform convolution in Fourier space
  y = ifft2(fft2(x) .* fft2(f, sx(1), sx(2)));
  y = cnv2slice(y, sf(1):sx(1), sf(2):sx(2));
elseif all(sx <= sf)  % x is smaller or equal than f
  y = cnv2(f, x);
else
  % x and f are incomparable
  error('[cnv2.m] x must be at least as large as f or vice versa.');
end
return

%%%%%%%%%%%%%%%%%
function f = cnv2tp(x, y)
sx = size(x);
sy = size(y);
% perform the linear convolution in Fourier space
if all(sx >= sy)
  sf = sx - sy + 1;
  f = cnv2slice(ifft2(conj(fft2(x)).*fft2(cnv2pad(y, sf))), 1:sf(1), 1:sf(2));
elseif all(sx <= sy)
  sf = sy + sx - 1;
  f = ifft2(conj(fft2(x, sf(1), sf(2))).*fft2(cnv2pad(y, sx), sf(1), sf(2)));
else  % x and y are incomparable
  error('[cnv2.m] x must be at least as large as y or vice versa.');
end
f = real(f);  
return

%%%%%%%%%%%%%
function B = cnv2pad(A, sf);
% PAD with zeros from the top-left
i = sf(1);  j = sf(2);
[rA, cA] = size(A);
B = zeros(rA+i-1, cA+j-1);
B(i:end, j:end) = A;
return

%%%%% END OF mbd.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
