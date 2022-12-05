function  pos_elements=generate_det_pos_file3D(NDet,AngCov,Radius,par,varargin)
% GENERATE_DET_POS_FILE3D Generates xml with detector positions, to be loaded
% with ViewMSOT
%
% imp_resp = generate_det_pos_file(NDet,AngCov,StepAngle,StartAngle,Radius,par,varargin)
% Generate position coordinate 
%
% Parameters:
% - NDet: Number of detector elements
% - AngCov: Angular coverage in degrees
% - Radius: Radius in m 
% - par: cup specific transducer parameters
% - filename: (optional)filename


% calcualted paramters
ang_hole=asin(par.R_in/(2*Radius)); %half opening angle of hole
ang_cov00=(0.5*AngCov/360)*2*pi; %angular coverage (half)
ang_cov0=ang_cov00-ang_hole; % half angular coverage minues hole
tot_size0=ang_cov0*Radius; %length of the detector curvature along elevation
row_num0=tot_size0/(par.el_size_elevation+par.gap); %number of rows before rounding 
row_num=round(row_num0); %number of rows after rounding
tot_size=(par.el_size_elevation+par.gap)*row_num; %length of the detector curvature along elevation for int number of elements
ang_cov=tot_size/Radius; %angular coverage minues hole for int number of elements
delta_theta=ang_cov/row_num; %step angle in elevation
theta0=pi/2+ang_hole+0.5*(par.el_size_elevation+par.gap)/Radius; %start angle in elevation (top) + offset due to hole
phi0=0; %start angle for azimuth 

% int variables
THETA_STEPS0=[theta0:delta_theta:theta0+(row_num-1)*delta_theta]; %vector with different angles in elevation
el_num_per_row0=zeros(size(THETA_STEPS0)); %number of elements per row
el_num_per_row=zeros(size(THETA_STEPS0)); %rounded number of elements per row 
el_azimuth_size_per_row=zeros(size(THETA_STEPS0)); %element size in azimuth after rounding
el_size_azimuth_vec=[]; %vector with element size in azimuth after rounding
THETA=[]; %vector with final elevation angles
PHI=[]; %vector with final azimuth angles

% loop over rows
for i=1:row_num
    theta_row=THETA_STEPS0(i); %angle of row
    Rcirc=Radius*sin(theta_row-pi/2); %radius of circle corresponding to a row
    el_num_per_row0(i)=2*pi*Rcirc/(par.el_size_azimuth+par.gap); %number of elements per row
    el_num_per_row(i)=round(el_num_per_row0(i)); %rounded number of elements per row 
    el_azimuth_size_per_row(i)=(2*pi*Rcirc)/el_num_per_row(i)-par.gap; %element size in azimuth after rounding 
    PHI0=[1:el_num_per_row(i)]*(2*pi/el_num_per_row(i)); %azimuth angles per row
    THETA0=ones(size(PHI0))*theta_row; %elevation angles per row
    PHI=cat(2,PHI,PHI0); %vector with final azimuth angles
    THETA=cat(2,THETA,THETA0); %vector with final elevation angles
    el_size_azimuth_vec=cat(2,el_size_azimuth_vec,el_azimuth_size_per_row(i)*ones(size(PHI0))); %vector with element size in azimuth after rounding
end

if numel(PHI) ~= NDet
   h1 = msgbox(['Calculated detector number (' num2str(numel(PHI)) ') does not match specification (' num2str(NDet) ')']);
   pause(2);
   delete(h1); 
   
   if numel(PHI) > NDet
      h1 = msgbox(['Increase el_size_azimuth']);
      pause(2);
      delete(h1);  
   else
      h1 = msgbox(['Decrease el_size_azimuth']);
      pause(2);
      delete(h1);  
   end
   return
end

% step angles
delta_phi_el=((el_size_azimuth_vec)/(Radius)); % vector with step angles in azimuth
delta_theta_el=((par.el_size_elevation)/(2*pi*Radius))*2*pi; % step angles in elevation

% transform to cartesian coordinates
pos_elements=zeros(length(PHI),3);
det_num_tot=length(PHI);
[pos_elements(:,1),pos_elements(:,2),pos_elements(:,3)] = sph2cart(PHI, THETA, ones(size(PHI))*Radius);


% plot
figure
plot3(0,0,0,'r.')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis image
hold
for j=1:NDet
plot3(pos_elements(j,1),pos_elements(j,2),pos_elements(j,3),'b*')
end
hold
axis equal
grid on

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
    

      
end


