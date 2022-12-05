classdef (Abstract) msotROIMeasure
    properties (Abstract, Constant)
        Description;        % Long Description
        Short;              % Short Description
        Ident;              % ID for adressing
        Priority;
    end
    
    methods (Abstract)
        mx = calculate(obj);
    end
    
    methods (Static)
        function obj = create(str)
%             if isempty
            
            
            
            try
                obj = eval(['msotROIMeasure' str '()']);
            catch ex
                if strcmp(ex.identifier,'MATLAB:UndefinedFunction')
                    error('Undefined ROI Measure (%s)',str);
                else
                    rethrow(ex);
                end
            end
        end
        
        function [c,s,d,i] = findDescendants()
           % find directory of class files
           dname = fileparts(which('msotROIMeasure')); 
           
           % iterate through all .m files
           flist = dir([dname filesep 'msotROIMeasure*.m']);
           s = {}; d = {}; i = {}; c = {}; p = []; cc = 0;
           for jj = 1:numel(flist),
               fn = strrep(flist(jj).name,'.m','');
               mc = meta.class.fromName(fn);
               % if superclass matches, 
               if ~isempty(mc.SuperclassList) && strcmp(mc.SuperclassList(1).Name,'msotROIMeasure'),
                   cc = cc + 1;
                   s{cc} = mc.PropertyList(ismember({mc.PropertyList.Name},'Short')).DefaultValue;
                   d{cc} = mc.PropertyList(ismember({mc.PropertyList.Name},'Description')).DefaultValue;
                   i{cc} = mc.PropertyList(ismember({mc.PropertyList.Name},'Ident')).DefaultValue;
                   p(cc) = mc.PropertyList(ismember({mc.PropertyList.Name},'Priority')).DefaultValue;
                   c{cc} = fn;
               end
           end
           
           [~,ind] = sort(p,'descend');
           c = c(ind);
           s = s(ind);
           d = d(ind);
           i = i(ind);
        end
    end
end



    
    