classdef msotLogger < handle
    % msotLogger Standard Logging Class
    %   Handles logging to central logfile and console. Specify an
    %   application Name to determine the log folder.
    %
    %   Per default, the log file is contained in the APPDATA\iThera
    %   folder, where a subfolder is created for each Application. Log File
    %   Names have the format <appName>_<date>-<time>.log, e.g.
    %       cLabs_20170104-180305.log
    %
    %   All messages that are lower than the specified limit are logged,
    %   these can be controlled separately for the different facilities.
    %
    %   In one App, only initialise once in a global variable, e.g.:
    %       >> global msotLOG;
    %       >> msotLOG = msotLogger(appName);
    %   Then, in subfunctions with separate namspaces or objects, use the
    %   global keyword again to access the global object:
    %       >> global msotLOG;
    %       >> msotLOG.info(myclass,'My message');
    %   
    %   NOTE: All log functions support sprintf style arguments
    %   
    %   Loglevels are:
    %   0 = error
    %   1 = warning
    %   2 = info
    %   3 = debug
    properties (SetAccess = private)
        FID@double scalar = double(0);  
        creation@datetime scalar;
        filename@char vector = char([]);
    end
    
    properties
        loglevel@int8 scalar = int8(2);           % Logging to file,   -1 = noOutput, 0 = error, 1 = warning, 2 = info, 3 = debug
        outlevel@int8 scalar = int8(3);           % Logging to STDOUT, -1 = noOutput, 0 = error, 1 = warning, 2 = info, 3 = debug
        appname@char vector = char([]);           % Name of Application to determine FileName
    end
    
    properties (Dependent)
        ID@uint16;
    end
    
    properties (Constant)
        sourceLength@uint8 = uint8(20);             % Length of Source Identifier
    end
    
    methods
%% Constructor / Destructor
        function obj = msotLogger(varargin)
            obj.creation = datetime;

            switch numel(varargin) 
                case 0, obj.appname = 'unknown'; obj.loglevel = int8(-1); obj.outlevel = int8(127);
                case 1, obj.appname = varargin{1}; obj.loglevel = int8(2); obj.outlevel = int8(3); 
                case 2, obj.appname = varargin{1}; obj.loglevel = varargin{2}; obj.outlevel = int8(3); 
                case 3, obj.appname = varargin{1}; obj.loglevel = varargin{2}; obj.outlevel = varargin{3}; 
                otherwise error('msotLogger:invalidConstructor','Invalid number of arguments to msotLogger Constructor (%i)',numel(varargin));
            end
            
            % If Logging to file is enabled, open the respective file
            if obj.loglevel >= 0
                obj.openLogFile();
            end
            
            obj.info(obj,'Logger initialised. Output Level = %i',obj.outlevel);
        end
        
        function delete(obj)
            obj.info(obj,'Terminating msotLogger...');
            obj.closeLogFile();
        end
        
%% Accessors
        % Unique ID to mimic msotHandle
        function id = get.ID(obj)
            id = GetMD5(datestr(obj.creation,'yymmddHHMMSSFFF'),'8bit','base64');
        end
        
        % ** Log Level modifications - tracked in Log...
        function set.loglevel(obj,ll)
            obj.info(obj,'LogLevel is now %i (was %i)',ll,obj.loglevel);
            obj.loglevel = ll;
            if ll >= 0 && obj.FID <= 0 %#ok<MCSUP>
                obj.openLogFile();
            end
        end

        % ** STDOUT Output Level modifications - tracked in Log...
        function set.outlevel(obj,ol)
            obj.info(obj,'OutLevel is now %i (was %i)',ol,obj.outlevel);
            obj.outlevel = ol;
        end
        
%% Log File 
        function fid = openLogFile(obj,varargin)
            % if a logfile is already open, close it
            isopen = true;
            try ftell(obj.FID); catch, isopen = false; end
            if obj.FID ~= 0 && isopen
                obj.info(obj,'Closing old LogFile (%i)',obj.FID);
                fclose(obj.FID);
                obj.FID = 0;
            end
            
            % determine file name from arguments or specified filename
            if numel(varargin) >= 1
                fname = varargin{1};
            elseif ~isempty(obj.filename)
                fname = obj.filename;
            % otherwise generate generic filename
            else
                dirname = [getenv('APPDATA') filesep 'iThera' filesep obj.appname filesep 'logs'];
                % create directory if it doesnt exist
                if ~exist(dirname,'dir')
                    res = mkdir(dirname);
                    if ~res, warning('msotLogger:invalidDir','Cannot create Output Directory %s',dirname); return; end
                end

                fname = [obj.appname '_' datestr(obj.creation,'yyyymmdd-HHMMSS') '.clog'];
                fname = [dirname filesep fname];
            end

            % open logfile
            [fid,msg] = fopen(fname,'w');
            if fid < 0
                obj.warn(obj,'msotLogger:cannotCreateFile','Cannot create file %s - %s.\nNot logging anything to file',fname,msg);
            else
                obj.info(obj,'Opened LogFile %s (Loglevel %i - Handle %i)',fname,obj.loglevel,fid);
                obj.filename = fname;
                obj.FID = fid;
            end
        end            
        
        function closeLogFile(obj)
            if obj.FID > 0
                try
                    fclose(obj.FID);
                end
                obj.FID = 0;
            end
        end
        
%% Actual Logging Functions
        function log(obj,source,level,msg,varargin)
            if isempty(obj)
                global msotLOG;  %#ok<TLEV>
                if isempty(msotLOG), msotLOG = msotLogger(); end 
                obj = msotLOG;
            end
            if level > obj.loglevel && level > obj.outlevel, return; end
            lev = '';
            switch level
                case 0, lev = 'error';
                case 1, lev = 'warn ';
                case 2, lev = 'info ';
                case 3, lev = 'debug';
                otherwise, warning('Invalid log level (%i)',level);
            end
            
            % source
            mc = metaclass(source);
            src = sprintf(['%-' num2str(obj.sourceLength) 's'],mc.Name);
            src = src(1:obj.sourceLength);
            
            % append ID
            if any(ismember({mc.PropertyList.Name},'sID'))
                id = source.sID;
            else
                id = uint16(0);
            end
            % For double and single use first four chars of base64
            if isa(id,'double') || isa(id,'single'), id = char(org.apache.commons.codec.binary.Base64.encodeBase64(typecast(id,'int8'))); id = id(1:4);
            % for Integers use HEX
            elseif isinteger(id), id = dec2hex(uint16(id),4); 
            % for chars just use 4 chars to make all the same length
            else, id = sprintf('%4s',id); id = id(1:4); 
            end
            
            if ~isempty(varargin) && isa(varargin{end},'MException')
                ex = varargin{end};
                varargin = varargin(1:end-1);
            else
                ex = MException.empty();
            end
            
            % format message
            msg = sprintf(msg,varargin{:});
            msg = strrep(msg,'\','\\');     % escape backslashes in file names
            msg = sprintf('%s ** %s ** [%s %s] %s\n',datestr(datetime,'yyyy-mm-dd HH:MM:SS'),lev,src,id,msg);
            
            % print to log and stdout
            if level <= obj.loglevel && obj.FID > 0, fprintf(obj.FID,msg); end
            if level <= obj.outlevel, fprintf(msg); end
            
            if ~isempty(ex)
                obj.log(source,level,'############################################################################');
                obj.log(source,level,'## Caught MATLAB Exception:');
                obj.log(source,level,'##   Message: %s',ex.message);
                obj.log(source,level,'##   Identifier: %s',ex.identifier);
                obj.log(source,level,'##   Stack:');
                for jj = 1:numel(ex.stack)
                obj.log(source,level,'##   - %s (%s, line %i)',ex.stack(jj).name,ex.stack(jj).file,ex.stack(jj).line);
                end    
                obj.log(source,level,'############################################################################');
            end
        end
        
        function error(obj,source,id,msg,varargin)
            % issue an error (will call error)
            obj.log(source,0,[msg ' [' id ']'],varargin{:});
            try
                msg = sprintf(msg,varargin{:}{:});
            catch
            end
            error(id,[id ': ' msg]);
        end
        
        function warn(obj,source,id,msg,varargin)
            % issue a warning (will call matlab warning)
            obj.log(source,1,[msg ' [' id ']'],varargin{:});
            try
                msg = sprintf(msg,varargin{:}{:});
            catch
            end
            warning(id,[id ': ' msg]);
        end
        
        function info(obj,source,msg,varargin)
            obj.log(source,2,msg,varargin{:});
        end
        
        function debug(obj,source,msg,varargin)
            obj.log(source,3,msg,varargin{:});
        end
        
    end
end