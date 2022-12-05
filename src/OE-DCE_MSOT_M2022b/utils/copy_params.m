function par = copy_params(par,cpar)
% par = copy_params(par,newpar)
%
% Copy parameters from struct "newpar" to struct par that is supplied with
% the defaults.

% transfer all fields to parameter array
fx = fieldnames(cpar);
for j = 1:numel(fx)
    par = setfield(par,fx{j},getfield(cpar,fx{j}));
end
