function  pos_elements=generate_det_pos_file3D_luis(NDet, AngCov, Radius, par, plot, varargin)
% GENERATE_DET_POS_FILE3D 
% Generates xml with detector positions, to be loaded with ViewMSOT
%
% imp_resp = generate_det_pos_file(NDet,AngCov,StepAngle,StartAngle,Radius,par,varargin)
% Generate position coordinate 
%
% Parameters:
% - NDet:   Number of detector elements
% - AngCov: Angular coverage in degrees
% - Radius: Radius in m 
% - par:    Cup specific transducer parameters
% - filename: (optional)filename


% Calcualted parameters
theta = zeros(NDet,1);
phi = zeros(NDet,1);
if numel(par.l_1) == 1 
    par.l_1 = repmat(par.l_1,[1,par.ring_num]);
end

l_1 = [];
k = 0;
for ring_count = 1:par.ring_num
    lp(ring_count) = par.l_0 + ((1+2*(ring_count-1))/2)*par.l_1(ring_count) + (ring_count-1)*par.gap;
    
    if ring_count>1
       if par.l_1(ring_count) ~= par.l_1(ring_count-1-k)
          k = k +1;
          
          if isempty(l_1) 
              l_1 = par.l_1(ring_count-1); 
              ring_last = ring_count-1; 
          end
          
          lp(ring_count) = par.l_0 + (((2*ring_last))/2)*par.l_1(ring_last) + (1+2*(mod(ring_count,ring_last)-1))/2*par.l_1(ring_count) + (ring_count-1)*par.gap; 
       end
    end
end

for ring_count = 1:par.ring_num
    theta(par.el_r{ring_count}) = pi/2-lp(ring_count)/Radius;
end

for ring_count = 1:par.ring_num
    pos_r = par.el_r{ring_count};
    Dphi = 2*pi/numel(pos_r);
        for ii = 1:numel(pos_r)
            phi(pos_r(ii)) = ((1+2*(ii-1))/2)*Dphi;
        end
end


pos_elements(:,1) = Radius.*cos(theta).*cos(phi);
pos_elements(:,2) = Radius.*cos(theta).*sin(phi);
pos_elements(:,3) = Radius.*sin(theta);


%% Plot

if plot == 1
    figure
    plot3(0,0,0,'r.')
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    axis image
    hold
    for j=1:NDet
    plot3(pos_elements(j,1),pos_elements(j,2),pos_elements(j,3),'b.')
    end
    hold
    axis equal
    grid on
end


%% Write positions to file

if ~isempty(varargin)    
    if isempty(varargin{1})
        filename = ['det_pos_NDet', num2str(NDet), '_AngCov', strrep(sprintf('%.1f',AngCov),'.','d'), 'deg_Radius', strrep(sprintf('%.1f',Radius/1e-3),'.','d'), 'mm'];
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

