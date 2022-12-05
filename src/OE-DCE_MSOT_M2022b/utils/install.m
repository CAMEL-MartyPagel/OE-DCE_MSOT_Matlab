pathStr = java.lang.String( path ) ;
os = java.lang.System.getProperty('os.name') ;
if( os.toLowerCase.startsWith( 'windows' ) )
    matlabPath = pathStr.substring( 0, pathStr.indexOf( ';' ) ) ;
else
    matlabPath = pathStr.substring( 0, pathStr.indexOf( ':' ) ) ;
end ;

cdStr = java.lang.String( cd ) ;
if( cdStr.endsWith( 'msotlib' ) )
    
    corelibDir = char( cdStr ) 
    
    % write startup.m into matlabPath folder 
    startup_file = java.io.File( matlabPath.concat( filesep ).concat( 'startup.m' ) ) 
    
    if( startup_file.exists() )
        bkp = java.io.File( matlabPath.concat( filesep ).concat( 'startup.m.bkp' ) ) ;
        startup_file.renameTo( bkp ) ;
    end ;
    
    fw = java.io.BufferedWriter( java.io.FileWriter( startup_file ) ) ;
    
    str1 = ['path( path, ''' corelibDir ''' ) ; '] ;
    fw.write( str1, 0, length( str1 ) ) ;
    fw.newLine() ;
    str2 = ['path( path, ''' corelibDir filesep 'recon' ''' ) ; '] ;
    fw.write( str2, 0, length( str2 ) ) ;
    fw.newLine() ;
    str2 = ['path( path, ''' corelibDir filesep 'utils' ''' ) ; '] ;
    fw.write( str2, 0, length( str2 ) ) ;
    fw.newLine() ;
%     str3 = ['cd ''' corelibDir ''' ; '] ;
%     fw.write( str3, 0, length( str3 ) ) ;
%     fw.newLine() ;
    str4 = 'setenv(''KMP_DUPLICATE_LIB_OK'', ''TRUE'')' ;
    fw.write( str4, 0, length( str4 ) ) ;
    fw.newLine() ;
    
    str5 =  [ 'javaaddpath ''' corelibDir filesep 'MSOTBeans' filesep 'xmlbeans-2.5.0' filesep 'lib' filesep 'xbean.jar'' ;' ] ;
    fw.write( str5, 0, length( str5 ) ) ;
    fw.newLine() ;
    
    str6 = [ 'javaaddpath ''' corelibDir filesep 'MSOTBeans' filesep 'patbeans.jar'' ;' ] ;
    fw.write( str6, 0, length( str6 ) ) ;
    fw.newLine() ;

    str6b = [ 'javaaddpath ''' corelibDir filesep 'MSOTBeans' filesep 'msotbeans.jar'' ;' ] ;
    fw.write( str6b, 0, length( str6b ) ) ;
    fw.newLine() ;

    
    str7 = [ 'javaaddpath ''' corelibDir filesep 'recon' filesep 'recon.jar'' ;' ] ;
    fw.write( str7, 0, length( str7 ) ) ;
    fw.newLine() ;
    
    fw.flush() ;
    fw.close() ;
    
    
    
end ;

clear all; 
% clc