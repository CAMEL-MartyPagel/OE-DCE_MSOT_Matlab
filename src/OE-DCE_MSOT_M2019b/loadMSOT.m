function datainfo  = loadMSOT( varargin )
% LOADMSOT  Load .msot META information
%   datainfo = LOADMSOT()       opens Open File Dialog to choose file
%   datainfo = LOADMSOT(file)   opens the specified file
% 
% This function loads the META information on the selected scan subfolder
% including information about both RECONs and MSPs performed on the scan
% raw data.
%
% Required Parameters:
%   <none>
%
% Optional:
%   1: filename    path and name of the file to load.
%
% Return Values:
%   1: datainfo    a structure containing the meta information, [] on error 
%
% Example:
%   datainfo = loadMSOT('Scan_1\Scan_1.msot');


%% parameters
if ( nargin == 0 )
    
    path = getappdata( 0, 'last_path' ) ;
    if( isempty( path ) )
        path = cd ;
    end ;
    fc = javax.swing.JFileChooser( java.io.File( path ) ) ;
    fc.setAcceptAllFileFilterUsed(0) ;
    fc.addChoosableFileFilter( javax.swing.filechooser.FileNameExtensionFilter(  'MSOT Header',  'msot' ) ) ;
    format = 'xml' ;
    
    if( fc.showOpenDialog( [] ) == javax.swing.JFileChooser.APPROVE_OPTION )
        path = char( fc.getSelectedFile().getParent() ) ;
        setappdata( 0, 'last_path', path ) ;
        filename = char( fc.getSelectedFile().getPath() ) ;
    else
        sigMat = [] ;
        datainfo = [] ;
        return ;
    end ;
elseif( nargin >= 1 )
    filename = varargin{1} ;
else
    datainfo = [];
    error( 'Not enough input arguments' ) ;
end ;

datainfo = msotData(filename);










