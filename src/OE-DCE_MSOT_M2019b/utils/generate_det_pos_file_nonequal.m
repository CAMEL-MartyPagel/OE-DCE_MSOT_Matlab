
function  pos_elements = generate_det_pos_file_nonequal(NDet1,NDet2,AngCov1,AngCov2,Radius1,Radius2,Pitch1,Pitch2,varargin)

NDet = 2*NDet2 + NDet1;
AngCov = 2*AngCov2 + AngCov1;
Pitch = (Pitch1 + Pitch2)/2;

%%%% Position of Transducer array element
d_th1 = Pitch1/Radius1;    % the correponding angle between two neighboring elements.
d_th2 = Pitch2/Radius2;    % the correponding angle between two neighboring elements.

max_theta = (NDet1-1)/2*d_th1; 

theta1 = (-1*max_theta + d_th1*[0:1:NDet1-1]); 
theta2 = (theta1(1) - d_th2*[NDet2-1:-1:0]); 
theta3 = (theta1(end) + d_th2*[0:1:NDet2-1]); 
theta = [ theta2  theta1  theta3  ];

chx_t = [Radius2*sin(theta2)  Radius1*sin(theta1)  Radius2*sin(theta3) ];
chz_t = [-Radius2*cos(theta2)  -Radius1*cos(theta1)  -Radius2*cos(theta3) ];

%% linear central part

chx_t(NDet2+1:NDet2+NDet1) = linspace(chx_t(NDet2+1),chx_t(NDet2+NDet1),NDet1);
chz_t(NDet2+1:NDet2+NDet1) = chz_t(NDet2)-0.0008;

%% detector positions in cartesian coordinates

pos_elements=zeros(NDet,3);
pos_elements(:,1) = chx_t;
pos_elements(:,2) = chz_t;

% flip outer segment
pos_elements(193:256,1) = pos_elements(256:-1:193,1);
pos_elements(193:256,2) = pos_elements(256:-1:193,2);

%
%pos_elements([1:64,193:256],1) = 100000;
%pos_elements([1:64,193:256],2) = -100000;

%pos_elements_new(1:128,1) = pos_elements(65:192,1);
%pos_elements_new(1:128,2) = pos_elements(65:192,2);
%pos_elements = pos_elements_new; NDet = 128;

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
    schz = strcat('Chz=0,1,',num2str(NDet),',1:',schz,num2str(pos_elements_us(NDet,2),digits),']');             
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
    
   
    
%     % open file
%     FID2 = fopen([filename,'.txt'],'w');
% 
%     % 
%     output_cell2 = {};
%     output_cell2{1} = strcat('Radius: ',num2str(Radius,'%.17f'));
%     output_cell2{2} = strcat('Step Angle: ',num2str(StepAngle,'%.17f'));
%     output_cell2{3} = strcat('Start Angle: ',num2str(StartAngle,'%.17f'));
%     
%     % write file
%     fprintf(FID2, '%s\n', output_cell2{:});
% 
%     % close file
%     fclose(FID2);

end