function datainfoarr  = listMSOT2( varargin )
% LISTMSOT  List .msot META information from Study Folder
%   datainfo = LISTMSOT()      List scans in current directory
%   datainfo = LISTMSOT(dir)   Fist scans in the specified directory
% 
% This function lists the META information on the selected study folder
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: dirname:    path of the study to load.
%   2: listing:     display listing (default true)
%
% Return Values:
%   1: datainfoarr: a structure containing the meta information, [] on error 
%
% Example:
%   datainfo = loadMSOT('Scan_1\Scan_1.msot');
%   datainfo = loadMSOT('Scan_1\Scan_1.msot', false);


%% parameters
verbose = true;
if ( nargin == 0 )
   studydir = '.'; 
elseif (nargin >= 2)
    studydir = varargin{1};
    verbose = varargin{2};
elseif( nargin >= 1 )
   studydir = varargin{1} ;
else
   datainfo = [];
   error( 'Not enough input arguments' ) ;
end ;



%% file list
cfiles = rdir([studydir '\**\' '*.mso*']);

wbar = waitbar(0,'Parsing MSOT Files...');
for ind = numel(cfiles):-1:1
    waitbar((numel(cfiles)-ind)/numel(cfiles),wbar);
    xml_filename = cfiles(ind).name;
    datainfoarr(ind) = msotData(xml_filename,false);
    datainfoarr(ind).progress = true;
end
close(wbar);


%%
if numel(datainfoarr) > 0,
    [tmp ind] = sort(arrayfun(@(x) x.CreationTime, datainfoarr));
    datainfoarr = datainfoarr(ind);
end

if verbose,
    for ind = 1:numel(datainfoarr)
        fprintf('%s\t\t%s\t%s\t%s\n',datainfoarr(ind).FolderName,datainfoarr(ind).CreationTimeTxt,datainfoarr(ind).TransducerType,datainfoarr(ind).Name);
    end
end

return


