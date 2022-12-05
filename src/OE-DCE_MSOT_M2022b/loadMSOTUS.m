function [US, ts, datainfo] = loadMSOTUS(varargin)
% LOADMSOTUS  Load MSOT Ultrasound Images
%   [US ts rep] = LOADMSOTUS()       
%           opens Open File Dialog to choose file, then loads US images
%   [US ts] = LOADMSOTUS(file)  
%           opens the specified file, then loads the US images
%   [US ts] = LOADMSOTUS(datainfo)   
%           opens the specified file associated with struct, then loads the US images
% 
% This function loads the US images from the selected scan
%
%
% Return values:
% US         An array containing the US images in 8D
%            1) x
%            2) y
%            3) frame number
% ts         timestamps in s from start

%% parameters
US = [];
ts = [];
datainfo = [];
selMat = [];

% if no arguments given
if (nargin == 0 || isempty(varargin{1}))
    datainfo = loadMSOT();
else
    datainfo = varargin{1};
end
if ~isstruct(datainfo) && ischar(datainfo)
    datainfo = loadMSOT(datainfo);
end

if numel(varargin) >= 2,
   selMat = varargin{2}; 
end

if ~datainfo.USpresent,
    warning('No US data present - exiting...');
    return;
end;

% if no selMat is given, use all
if isempty(selMat),
    selMat = 1:numel(datainfo.US.timestamps);
else 
    % selMat starts at 0, but we want to have it starting at 1!
    selMat = selMat + 1;
end;

%% Initialise Progress Bar
wbar = waitbar(0,'Initialising, please wait...' ) ;


%% image return vector
n = double(datainfo.US.n);
if datainfo.US.nz ~= 0,  % for v4 Field of view definition
    nz = double(datainfo.US.nz);
else  % legacy case
    nz = n;
end
ts = datainfo.US.timestamps(selMat);
US = nan([nz n numel(selMat)],'single');
selMat = double(selMat);

%% Load Images
fname = strrep(datainfo.FileName,'.bin','.us');
[FID str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary US File ' fname ': ' str]);
    return;
end
% handle progress bar
waitbar(0,wbar,'Loading Reconstructions...');

% load all selected images sequentially
selMat = selMat(:);
for jj = 1:numel(selMat)
    id = selMat(jj);
    waitbar(jj/numel(selMat),wbar);
    
    ts(id) = datainfo.US.timestamps(id);
    try
        startbyte = (id-1)*n*nz*4;    %determine start position (single prec)
        fseek(FID,startbyte,-1);     % move towards position
        US(:,:,jj) = fread(FID,[nz n],'single');
    catch ex
        warning(['Cannot Read ReconFrame ' num2str(id) ', skipping...']);
    end
   
end
fclose(FID);

close(wbar) ;
