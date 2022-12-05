function [spec wl] = loadXVUEspectrum(fnames,varargin)

wls = [];
if numel(varargin) >= 1
    wls = varargin{1};
end

if ~iscell(fnames)
    fnames = {fnames};
end

wls = reshape(wls,1,numel(wls));
arr = msotSpectrum.create(fnames,wls);
spec = cell2mat(cellfun(@(x) x(:),arr,'UniformOutput',false));

% for c = 1:numel(fnames)
%     fname = [ 'C:\ProgramData\iThera\MSOTSystem\Spectra\' fnames{c} '.csv'];
%     if (~exist(fname,'file')),
%         fname = [ 'C:\ProgramData\iThera\MSOTSystem\Factory Spectra\' fnames{c} '.csv'];
%     end
% 
%     if (~exist(fname,'file')),
%         error('Cannot find Spectra of that name - filename: %s',fname);
%     end
%     
%     sfile = importdata(fname,';');
%     
%     clear wl abs;
%     for l = 1:numel(sfile)
%         tmp = textscan(strrep(sfile{l},',','.'),'%f;%f');
%         wl(l) = tmp{1};
%         abs(l) = tmp{2};
%     end
%     
%     if isempty(wls), 
%         spec(:,c) = abs; 
%     else 
%         spec(:,c) = spline(wl,abs,wls);
%         wl = wls;
%     end
% end