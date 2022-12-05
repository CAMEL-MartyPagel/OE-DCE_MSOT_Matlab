function [R, zpos, ts, datainfo] = loadMSOTMsp(varargin)
% LOADMSOTMsp  Load MSOT MSPs
%   [R wls zpos ts datainfo] = LOADMSOTMSP()       
%           opens Open File Dialog to choose file, then loads the first
%           Msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(file)  
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo)   
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum)   
%           opens the specified file, then loads the specified msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum,selMat)   
%           opens the specified file, then loads the specified msp
%           restricted to selMat, which is a multidimensional structure of 
%           frame ids. For all, use datainfo.MspNode(i).Structure
% 
% This function loads the MSOT MSPs from the selected scan folder. 
%
%
% Return values:
% R          An array containing the MSPed images in 8D
%            1) x
%            2) y
%            3) z (only if 3D detector)
%            4) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            5) z-position (if 2D system with translation stage)
%            6) REPETITION (currently unused)
%            7) Wavelength (see wls vector)
%            8) Individual Images (only if not averaged)

% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same 
%            dimension as R

%% parameters
mNode = 0;
selMat = [];
datainfo = [];

% if no arguments given
if (nargin == 0 || isempty(varargin{1}))
    datainfo = loadMSOT();
elseif isa(varargin{1}, 'msotData')
    datainfo = varargin{1};
else 
    datainfo = loadMSOT(varargin{1});
end

if isempty(datainfo) || ( ~isstruct(datainfo) && ~isa(datainfo,'msotData') )
    datainfo = loadMSOT(datainfo);
    %return
end

% if no reconNode number is given, use the first
if nargin >= 2
    mNode = varargin{2};
else
    mNode = 1;
end

if (numel(datainfo.MSPNode) < mNode)
    error('MSP not found');
end

% if no selMat is given, use all
if nargin >= 3
    selMat = varargin{3};
else
    selMat = datainfo.MSPNode(mNode).Structure;
end

%% Initialise Progress Bar
wbar = waitbar(0,'Initialising, please wait...' ) ;


%% initialise return array
sel = ~isnan(selMat);   % this will retain structure, marking only not-nans
svec = size(sel);       % retain structural information
selMat = selMat(sel);   % remove nans and make selMat 1D

% image return vector
nx = double(datainfo.MSPNode(mNode).Resolution);
ny = double(datainfo.MSPNode(mNode).ResolutionY);
if ny == 0  % case without FIELD-OF-VIEW, so pre 4.0
    ny = nx;
end

if datainfo.is3D
    nz = double(datainfo.MSPNode(mNode).ZResolution);
end

% initial vectors for zpositions and timestamps
zpos = zeros(size(selMat));
ts = zeros(prod(svec),1);

%% Load Images
fname = [datainfo.RealPath filesep 'MSPs' filesep datainfo.MSPNode(mNode).GUID '.bin'];
[FID, str] = fopen(fname,'r','n');
if (FID == -1)
    error(['Error opening binary MSP File ' fname ': ' str]);
    return
end

% Adapted from "PS 2020-02-26": detect precision of binary file (double or single)
% DataModelVersion added for objects of type MSPNodes
% if DataModelVersion exist and >=1.0, fileprecision='single'
if ~ischar(datainfo.MSPNode(mNode).DataModelVersion) || ~isprop(datainfo.MSPNode(mNode), 'DataModelVersion')
    
    fseek(FID, 0, 'eof');
    filesize = ftell(FID);

    selMat_tot = datainfo.MSPNode(mNode).Structure;
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
    dataModelVer = str2double(datainfo.MSPNode(mNode).DataModelVersion);
    
    if dataModelVer >= 1.0
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
waitbar(0,wbar,'Loading MSP...');

% load all selected images sequentially
for j = 1:numel(selMat)
    waitbar(j/numel(selMat),wbar);
    id = double(selMat(j));
    if id == 0 
        warning('Skipping frame...'); 
        continue 
    end
    slid = datainfo.MSPNode(mNode).SliceIndex(id);
    
    % copy Meta-Information from Scan Frame
    zpos(j) = datainfo.MSPNode(mNode).Slices(slid).ZPos;
    ts(j) = datainfo.MSPNode(mNode).Slices(slid).RelTime;
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
            if isprop(datainfo.MSPNode(mNode), 'DataModelVersion')
                if ~isempty(datainfo.MSPNode(mNode).DataModelVersion)
                    dataModelVer = str2double(datainfo.MSPNode(mNode).DataModelVersion);
                        if dataModelVer >= 1.0
                            R(:,:,j) = fliplr(transpose(squeeze(R(:, :, j))));
                        end
                end
            end
        end
        
    catch ex
        warning(['Cannot Read MspFrame ' num2str(id) ', skipping...']);
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
zpos = unique(zpos);

close(wbar) ;
