classdef (Abstract) msotHandle < handle
    properties (SetAccess = private)
        creation@datetime scalar = datetime.empty;
        ID@char vector = '';
        logger@msotLogger = msotLogger.empty();
    end
    
    properties (Dependent)
        sID@char;                   % short version of unqiue ID (3byte, 4char)
    end
    
    methods
        function obj = msotHandle()
            obj.creation = datetime;
            obj.ID = GetMD5([dec2hex(uint16(rand*2^16),4) datestr(obj.creation,'yymmddHHMMSSFFF')],'8bit','base64');
            global msotLOG;
            if isempty(msotLOG)
                msotLOG = msotLogger();
            end
            obj.logger = msotLOG;
            msotLOG.debug(obj,'Instance Created...');
        end
        
        function delete(obj)
            obj.logger.debug(obj,'Destroyed...');
        end
        
        function sid = get.sID(obj)
            sid = obj.ID;
            sid = sid(1:4);
        end
        
        function hash = generateHashCode(obj)
            % Generate Hash Code based on char conversion of obj. Should be
            % stable also across versions of the object, since properties 
            % set to default are not contained in the char representation
            % by
            hash = GetMD5(char(obj,true),'8bit','base64');
            hash = hash(1:8);
        end
    end
end