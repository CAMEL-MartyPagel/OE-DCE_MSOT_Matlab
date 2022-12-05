function  pos_elements=generate_det_pos_file(NDet,AngCov,Radius,varargin)
% GENERATE_DET_POS_FILE Generates xml with detector positions, to be loaded
% with ViewMSOT
%
% imp_resp = generate_det_pos_file(NDet,AngCov,StepAngle,StartAngle,Radius,varargin)
% Generate position coordinate 
%
% Parameters:
% - NDet: Number of detector elements
% - AngCov: Angular coverage in degrees
% - Radius: Radius in m 
% - filename: (optional)filename
% - clinical: true(default)


%% parameters
if ~isempty(varargin) && numel(varargin)==2
   clinical = varargin{2};
else clinical = 1; % default
end

StepAngle = ( AngCov*(pi/180) )/NDet; % step angle

if clinical
    StartAngle =  - 0.5*AngCov*(pi/180) + StepAngle/2  - pi/2; % start angle
else
    StartAngle =  - 0.5*AngCov*(pi/180) + StepAngle/2  - pi/2 + pi; % start angle
end
%% detector positions in spherical coordinates

pos_elements_sph=zeros(NDet,3);
for i=1:NDet
    pos_elements_sph(i,:)=[0,StartAngle+(i-1)*StepAngle,Radius];
end

%% detector positions in cartesian coordinates

pos_elements=zeros(NDet,3);
[pos_elements(:,1),pos_elements(:,3),pos_elements(:,2)]=sph2cart(pos_elements_sph(:,1),pos_elements_sph(:,2),pos_elements_sph(:,3));

% ensure symmetry
pos_elements_new = pos_elements;
pos_elements_new(NDet/2+1:NDet,2) = flipud((pos_elements_new(1:NDet/2,2)));
pos_elements_new(NDet/2+1:NDet,1) = flipud(-pos_elements_new(1:NDet/2,1));
pos_elements = pos_elements_new;


% print in .msot format
% digits = 17;
% for ii=1:NDet
%     sprintf(['<PROJECTION wavelength-ref="all" number="' num2str(ii) '"> \n'...
%         '<VALUE axis-ref="1">' num2str(pos_elements(ii,1),digits) '</VALUE> \n'...
%         '<VALUE axis-ref="2">' num2str(pos_elements(ii,2),digits) '</VALUE> \n'...
%         '<VALUE axis-ref="3">0</VALUE> \n'...
%         '</PROJECTION>'])
% end


%% display for ultrasound .ini file
pos_elements_us = pos_elements;
digits = 17;
schx='[';
schz='[';
for ii=1:NDet-1
    schx = strcat(schx,num2str(pos_elements_us(ii,1),digits),',');
    schz = strcat(schz,num2str(-pos_elements_us(ii,2),digits),',');  
end
    schx = strcat('Chx=0,1,',num2str(NDet),',1:',schx,num2str(pos_elements_us(NDet,1),digits),']');
    schz = strcat('Chz=0,1,',num2str(NDet),',1:',schz,num2str(-pos_elements_us(NDet,2),digits),']');             
% copy to [RUCT] field of ini file
disp(schx);
disp(schz);

%% write positions to file

if isempty(varargin)
    
else    
    if isempty(varargin{1}),
        filename = ['det_pos_NDet', num2str(NDet),'_AngCov',strrep(sprintf('%.1f',AngCov),'.','d'),'deg_Radius',strrep(sprintf('%.1f',Radius/1e-3),'.','d'),'mm'];
    else
        filename = varargin{1};
    end

    % open file
    FID = fopen([filename,'.xml'],'w');

    % generate header
    output_cell = {};
    output_cell{1}='<?xml version="1.0" encoding="UTF-8" standalone="yes"?>';
    output_cell{2}='<DataModelSensorPositions xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">';

    % generate detector positions 
    output_cell{3}='	<Sensors>';
    for i=1:size(pos_elements,1)

        pos_cell{1}='		<SensorPoint>';
        pos_cell{2}=strcat('			<x>',num2str(pos_elements(i,1),'%.17f'),'</x>');
        pos_cell{3}=strcat('			<y>',num2str(pos_elements(i,2),'%.17f'),'</y>');
        pos_cell{4}=strcat('			<z>',num2str(pos_elements(i,3),'%.17f'),'</z>');
        pos_cell{5}='		</SensorPoint>';
        output_cell=[output_cell,pos_cell];
        clear pos_cell
    end
    output_cell=[output_cell,'    </Sensors>'];
    output_cell=[output_cell,'</DataModelSensorPositions>'];

    % write file
    fprintf(FID, '%s\n', output_cell{:});

    % close file
    fclose(FID);
    
   
    
    % open file
    FID2 = fopen([filename,'.txt'],'w');

    % 
    output_cell2 = {};
    output_cell2{1} = strcat('Radius: ',num2str(Radius,'%.17f'));
    output_cell2{2} = strcat('Step Angle: ',num2str(StepAngle,'%.17f'));
    output_cell2{3} = strcat('Start Angle: ',num2str(StartAngle,'%.17f'));
    
    % write file
    fprintf(FID2, '%s\n', output_cell2{:});

    % close file
    fclose(FID2);
      
end


