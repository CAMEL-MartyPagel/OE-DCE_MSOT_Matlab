function [R, wls, zpos, ts, datainfo] = loadMSOTRecon(varargin)
% LOADMSOTRECON  Load MSOT Reconstructions
%   [R wls zpos ts datainfo] = LOADMSOTRECON()       
%           opens Open File Dialog to choose file, then loads the first
%           Recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(file)  
%           opens the specified file, then loads the first recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo)   
%           opens the specified file, then loads the first recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo,reconNum)   
%           opens the specified file, then loads the specified recon
%   [R wls zpos ts datainfo] = LOADMSOTRECON(datainfo,reconNum,selMat)   
%           opens the specified file, then loads the specified recon
%           restricted to selMat, which is a multidimensional structure of 
%           frame ids. For all, use datainfo.ReconNode(i).ReconStructure
% 
% This function loads the MSOT Reconstructions from the selected scan
% folder. 
%
%
% Return values:
%%% 2D data %%%
% R          An array containing the reconstructed images in 8D
%            1) x
%            2) y
%            3) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            4) z-position (if 2D system with translation stage)
%            5) REPETITION (currently unused)
%            6) Wavelength (see wls vector)
%            7) Individual Images (only if not averaged)

%%% 3D data %%%
% R          An array containing the reconstructed images in 8D
%            1) x
%            2) y
%			 3) z
%            4) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            5) z-position (if 2D system with translation stage)
%            6) REPETITION (currently unused)
%            7) Wavelength (see wls vector)
%            8) Individual Images (only if not averaged)


% wls        Wavelength Vector
% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same 
%            dimension as R

%%%% File history %%%%
% PS 2020-02-26: detect precision of binary file (double or single)

%% parameters
rNode = 0;
selMat = [];
datainfo = [];

% if no arguments given
if (nargin == 0 || isempty(varargin{1}))
    datainfo = loadMSOT();
else
    datainfo = varargin{1};
end
if ~isstruct(datainfo) && ischar(datainfo)
    datainfo = loadMSOT(datainfo);
end

% if no reconNode index is given, use the first
if nargin >= 2
    rNode = varargin{2};
else
    rNode = 1;
end

if (numel(datainfo.ReconNode) < rNode)
    error('Reconstruction not found');
end

% if no selMat is given, use all
if nargin >= 3 && ~isempty(varargin{3})
    selMat = varargin{3};
else
    selMat = datainfo.ReconNode(rNode).ReconStructure;
end

%% Initialise Progress Bar
wbar = waitbar(0,'Initialising, please wait...' ) ;


%% initialise return array
sel = ~isnan(selMat);   % this will retain structure, marking only not-nans
svec = size(sel);       % retain structural information
selMat = selMat(sel);   % remove nans and make selMat 1D

% image return vector
nx = double(datainfo.ReconNode(rNode).Resolution);
ny = double(datainfo.ReconNode(rNode).ResolutionY);
if ny == 0  % case without FIELD-OF-VIEW, so pre 4.0
    ny = nx;
end
% PS 2020-02-26: detect resolution in z if 3D data
if datainfo.is3D
    nz = double(datainfo.ReconNode(rNode).ZResolution);
end

% initial vectors for wavelengths and zpositions and timestamps
wls = zeros(size(selMat));
zpos = zeros(size(selMat));
ts = zeros(prod(svec),1);

%% Load Images
fname = [datainfo.RealPath filesep 'RECONs' filesep datainfo.ReconNode(rNode).GUID '.bin'];
[FID, str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary Recon File ' fname ': ' str]);
    return
end

% PS 2020-02-26: detect precision of binary file (double or single)
% if dataModelVer exist and >=2.0, fileprecision='single'
if ~ischar(datainfo.ReconNode(rNode).DataModelVersion) || ~isprop(datainfo.ReconNode(rNode), 'DataModelVersion')

    fseek(FID, 0, 'eof');
    filesize = ftell(FID);
    
    selMat_tot = datainfo.ReconNode(rNode).ReconStructure;
    sel_tot = ~isnan(selMat_tot);
    svec_tot = size(sel_tot); % To reliably estimate the fileprecision from the number of pixels encoded in the .bin, the total number of frames has to be considered for the calculation -- not only the frames loaded by the user!
    
    if datainfo.is3D
        total_pixelnum = nx*ny*nz*prod(svec_tot);
    else
        total_pixelnum = nx*ny*prod(svec_tot);
    end
    byte_pixelsize = filesize/total_pixelnum;
    if byte_pixelsize > 5  
        fileprecision = 'double';   % double=8
    else
        fileprecision = 'single';   % single=4
    end
    
else
    dataModelVer = str2double(datainfo.ReconNode(rNode).DataModelVersion);
    
    if dataModelVer >= 2.0
        fileprecision = 'single';
    else
        fileprecision = 'double';        
    end
    
end

if datainfo.is3D
    R = nan([nx ny nz prod(svec)],fileprecision);
else
    R = nan([nx ny prod(svec)],fileprecision);
end

% handle progress bar
waitbar(0,wbar,'Loading Reconstructions...');

% load all selected images sequentially
for j = 1:numel(selMat)
    waitbar(j/numel(selMat),wbar);
    id = double(selMat(j));
    if id == 0 || isnan(id)
        fprintf('Skipping frame %i - id is 0\n',j);
        continue
    end
    
    % copy Meta-Information from Scan Frame
    wls(j) = datainfo.ReconNode(rNode).Frames(id).Wavelength;
    zpos(j) = datainfo.ReconNode(rNode).Frames(id).ZPos;
    ts(j) = datainfo.ReconNode(rNode).Frames(id).RelTime;
    try
        if datainfo.is3D
            if strcmp(fileprecision, 'double')
                startbyte = (double(id)-1)*double(nx)*double(ny)*double(nz)*8;    %determine start position (double prec), one value is counted for each z slice
            else
                startbyte = (double(id)-1)*double(nx)*double(ny)*double(nz)*4;
            end
            
            fseek(FID,startbyte,-1);     % move towards position
            for zi = 1:nz
                R(:,:,zi,j) = fread(FID,[nx ny],fileprecision);
            end
            
        else
            if strcmp(fileprecision, 'double')
                startbyte = (double(id)-1)*double(ny)*double(nx)*8;    %determine start position (double prec)
            else
                startbyte = (double(id)-1)*double(ny)*double(nx)*4;
            end
            
            fseek(FID,startbyte,-1);     % move towards position
            R(:,:,j) = fread(FID,[nx ny],fileprecision);
            
            % Account for different axes orientation in latest SW versions
            if isprop(datainfo.ReconNode(rNode), 'DataModelVersion')
                if ~isempty(datainfo.ReconNode(rNode).DataModelVersion)
                    dataModelVer = str2double(datainfo.ReconNode(rNode).DataModelVersion);
                        if dataModelVer >= 2.0
                            R(:,:,j) = fliplr(transpose(squeeze(R(:, :, j))));
                        end
                end
            end
        end
        
    catch ex
        warning(['Cannot Read ReconFrame ' num2str(id) ', skipping...']);
    end
   
end
fclose(FID);

%% reshape data to fit requirements
if datainfo.is3D
    R = reshape(R,[nx ny nz svec]);
else
    R = reshape(R,[nx ny svec]);
end
ts = reshape(ts,svec);
wls = unique(wls);
zpos = unique(zpos);

close(wbar) ;
