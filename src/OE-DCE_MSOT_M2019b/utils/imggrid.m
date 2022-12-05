function ax = imggrid(varargin)
% IMGGRID   Create Grid of Axis 
% 
% Parameters:
% - Figure Handle (optional, otherwise current figure)
% - Number of Columns
% - Number of Rows
% - Size of Images (in px, 1 or 2 elements)
% - Spacing of Images (optional, default 0)
% - Number of Images (optional, default rows*cols)


ai = 1;
if nargin >= 1 && ishandle(varargin{ai}),
    fig = varargin{ai};
    ai = ai + 1;
else
    fig = gcf;
end

if nargin-ai+1 < 3, error('Not enough input parameters'); end

cols = varargin{ai}; ai = ai + 1;
rows = varargin{ai}; ai = ai + 1;
n = varargin{ai}; ai = ai + 1;

if nargin >= ai, spac = varargin{ai}; ai = ai +1; else spac = 0; end
if nargin >= ai, nimg = varargin{ai}; ai = ai +1; else nimg = rows*cols; end

if numel(n) == 1, n = [n n]; end

%% Adapt Window Position
wi = cols*n(1)+(cols+1)*spac;
he = rows*n(2)+(rows+1)*spac;
fig.Units = 'pixels';
fig.Position = [100 100 wi he];


%% Create Axis
ii = 0; col = 0; row = 1;
p(4) = n(2); p(3) = n(1);
while (ii < nimg),
    ii = ii + 1;
    col = col + 1;
    if col > cols, row = row + 1; col = 1; end
    
    p(1) = col*spac+(col-1)*n(1);
    p(2) = (rows-row+1)*spac+(rows-row)*n(2);
    ax(ii) = axes('Units','pixels','Position',p);
    axis image; axis off;
end
