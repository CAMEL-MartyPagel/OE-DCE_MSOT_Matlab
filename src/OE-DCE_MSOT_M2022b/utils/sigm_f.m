function [ map,s_shape ] = sigm_f( x,varargin )
slope=1; 
center=0.5;
max_val=1;

if (numel(varargin) >= 1),
    slope = varargin{1};
end
if (numel(varargin) >= 2)
    center= varargin{2};
end
if (numel(varargin) >= 3),
    max_val= varargin{3};
end


map=max_val./(1+exp(-slope.*(x-center)));

l=linspace(0,1,255);
s_shape=max_val./(1+exp(-slope.*(l-center)));


end

