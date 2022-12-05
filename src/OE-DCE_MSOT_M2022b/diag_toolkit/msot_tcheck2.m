%% input

% source type
%TCheck.source_type='ViewMSOT'
TCheck.source_type='LabTools2'

% data path
%TCheck.path='E:\MSOT Data\11408_1009&1010\Scan_1\Scan_1.msot';
%TCheck.path = 'E:\MSOT Data\Imasonic11408';
%TCheck.path = 'E:\MSOT Data\9967_1009&1010\'
TCheck.path = 'E:\MSOT Data\9967_1009&1010\170606_16-55-55_scan'

% simulation
%TCheck.sim_file= 'P:\_Alex\Code\com.ithera-medical.msotlib\transducer simulation\detinVision512_Rcurve37(2016-12-21-16-04-16)';
TCheck.sim_file='';

% probe file
%TCheck.probe_file = 'J:\Operations\Purchasing\Imasonic\1_UT Arrays 2D\11408 (InVision512)\Software Input Files\inVision512_Rcurve37.probe'
TCheck.probe_file = 'J:\Operations\Purchasing\Imasonic\1_UT Arrays 2D\09967 (inVision256-TF)\Software Input Files\InVision256-TF.probe'

% UT serial number
TCheck.UTSN='1009-1010';
%TCheck.UTDID='Imasonic 11408';
TCheck.UTDID='Imasonic 9967';
 
% channel offset
TCheck.Misc.ChannelOffset=0;

% sensor position modification
TCheck.Misc.RRoff=0;
TCheck.Misc.AZoff=pi;

% set filter parameters
TCheck.filter.f_LPF=10e6;
TCheck.filter.f_HPF=0.05e6;
TCheck.filter.fs=40e6;
TCheck.filter.order=4;
TCheck.filter.tukeywin_r=800;

% center of mass threshhold
TCheck.CenterOfMass.thresh=0.90;

% temperature
TCheck.Temp=22

% reconstruction
TCheck.Recon.fov_size = 35e-3; % size of the field of view (one axis)
TCheck.Recon.n = 2000; % pixel number for reconstruction (one axis)
TCheck.Recon.fmin = 0.05; % filter limits (MHz)
TCheck.Recon.fmax = 7;

%%
%load probe class
TCheck.probe=msotProbe.readXML(TCheck.probe_file);

%% load signals

if strcmp(TCheck.source_type,'LabTools2')
    sigMat=zeros(TCheck.probe.NumDetectors,2032,81*61,'uint16');
   
    for i=1:61
        i-1
        sigMat(:,:,1+(i-1)*81:81+(i-1)*81) = permute(readNPY([TCheck.path,'\data_slice_',num2str(i),'.npy']), [2 1 3]);
    end
           
    posx = linspace(-15,15,9);
    posy = linspace(-15,15,9);
    posz = linspace(-6,6,61);
    
    [X,Y,Z] = meshgrid(posx,posy,posz);
    
%     posx=X(:)';
%     posy=Y(:)';
%     posz=Z(:)';
    
    posx=Y(:)';
    posy=-X(:)';
    posz=Z(:)';
    
    TCheck.filter.skip_channel=3;
    
    info.zStep_scan = 61-1;
    info.xStep_scan = 9-1;
    info.yStep_scan = 9-1;
    
    info.zMin_scan=-6;
    info.zMax_scan=6;
    info.zCenter_scan=0;   
    
    sigMat = permute(sigMat, [2 1 3]);    
end

if strcmp(TCheck.source_type,'ViewMSOT')
    
    dinfo = loadMSOT([TCheck.path,'Scan_1\Scan_1.msot']);
    
    posx = linspace(-15,15,9);
    posy = linspace(-15,15,9);
    posz = linspace(-6,6,62);
    [X,Y,Z] = meshgrid(posx,posy,posz);
    posx=X(:)';
    posy=Y(:)';
    posz=Z(:)';
    
    par.usePower = 0;
    par.average = 1;
    sigMat = loadMSOTSignals(dinfo,[],par);
    sigMat=reshape(sigMat,2030,dinfo.HWDesc.NumDetectors,81*62);
    
    TCheck.filter.skip_channel=1;
    info.zStep_scan = 62-1;
    info.xStep_scan = 9-1;
    info.yStep_scan = 9-1;
    
    info.zMin_scan=-6;
    info.zMax_scan=6;
    info.zCenter_scan=0;   
end



%% load simulation

if ~isempty(TCheck.sim_file) %check if sim file is present
    
    load(TCheck.sim_file) % load file
    
    out.sig_amp_mat=flip(out.sig_amp_mat,4); % FLIP ORDER OF CHANNELS
       
    % grid vector for simulation
    poszz=(-(par.Nsteps_x-1)/2:+(par.Nsteps_x-1)/2)*par.stepx*1e3;
    posxx=(-(par.Nsteps_y-1)/2:+(par.Nsteps_y-1)/2)*par.step*1e3;
    posyy=(-(par.Nsteps_y-1)/2:+(par.Nsteps_y-1)/2)*par.step*1e3;

    % grid for interpolation (simulation)
    [XX,YY,ZZ] = meshgrid(posxx,posyy,poszz);
    
else
    
    out.sig_amp_mat=-1;
    
end
%%


%% detector positions

% convert to spherical
[AZ,EL,RR] = cart2sph(TCheck.probe.Sensors(:,1),TCheck.probe.Sensors(:,2),TCheck.probe.Sensors(:,3));

% change radius
RR=RR+TCheck.Misc.RRoff;

% change in-plane orientation
AZ=AZ+TCheck.Misc.AZoff; 

% convert back to cartesian
[sensorX,sensorY,sensorZ] = sph2cart(AZ,EL,RR);
TCheck.sensorPOS=[sensorX,sensorY,sensorZ];

%% filter signals and apply window

% obtain filter coefficients
[b,a] = butter(TCheck.filter.order,[TCheck.filter.f_HPF/(TCheck.filter.fs/2) TCheck.filter.f_LPF/(TCheck.filter.fs/2)],'bandpass');

% window
L = 2030;
TWin = tukeywin(L,TCheck.filter.tukeywin_r/L); 

% init
TCheck.sigMatFilt=zeros(size(sigMat(TCheck.filter.skip_channel:end,:,:)),'single');

% filter
for i=1:TCheck.probe.NumDetectors
    i
    for j=1:size(sigMat,3) 
        TCheck.sigMatFilt(:,i,j) = TWin.*filtfilt( b, a, -double(sigMat(TCheck.filter.skip_channel:end,i+TCheck.Misc.ChannelOffset,j) ));
    end
end


%% determine and select center microsphere in XY and XYZ

% sort distance of recorded XY position to center
[~,I]=sort(sqrt(posx.^2+posy.^2));

% select the signal matrices correponding to the center position in XY
TCheck.sigMatFilt_XYcenter=TCheck.sigMatFilt(:,:,I(1:info.zStep_scan+1));

% determine center in Z
meanAmpVecCenter=mean(squeeze(max(TCheck.sigMatFilt_XYcenter)-min(TCheck.sigMatFilt_XYcenter)));
maxSliceNum=find(meanAmpVecCenter==max(meanAmpVecCenter));

%select center signal matrix
TCheck.sigMatFilt_XYZcenter=TCheck.sigMatFilt_XYcenter(:,:,maxSliceNum);

% sort distance of recorded XY position to center
[~,I]=sort(abs(posz));

% select the signal matrices correponding to the center position in Z
TCheck.sigMatFilt_Zcenter=TCheck.sigMatFilt(:,:,I(1:2*(info.xStep_scan+1)*(info.yStep_scan+1)));


%% generate sensitivity field (3D matrix)
TCheck.AmpMat=zeros(size(X,1),size(X,2),size(X,3),TCheck.probe.NumDetectors);
TCheck.AmpVec=zeros(TCheck.probe.NumDetectors,size(TCheck.sigMatFilt,3));

for j=1:TCheck.probe.NumDetectors
    j
    T=squeeze(TCheck.sigMatFilt(:,j,:));
    TCheck.AmpVec(j,:)=max(T,[],1)-min(T,[],1);
    F = scatteredInterpolant(posx',posy',posz',double(TCheck.AmpVec(j,:)'));
    %F = scatteredInterpolant(posx,posy,posz,double(TCheck.AmpVec(j,:)'));
    TCheck.AmpMat(:,:,:,j)=F(X,Y,Z); 
end

TCheck.AmpMat=permute(TCheck.AmpMat, [2,1,3,4]); % FLIP X AND Y !!!!!!!!!!!!!

%% sensitifity field orientation

% intit
TCheck.MomentsOfInertia.IMat=zeros(3,3,TCheck.probe.NumDetectors);
if ~isempty(TCheck.sim_file) %check if sim file is present
    TCheck.MomentsOfInertia.IMat_sim=zeros(3,3,TCheck.probe.NumDetectors);
end

% generate output path
mkdir([TCheck.path,'\figures\orientation\'])


%moments of inertia and center of mass
for k=1:TCheck.probe.NumDetectors
    
    k
    
    % indiviual moments
    TM=TCheck.AmpVec(k,:);
    Ixx=sum(TM.*(posy.^2+posz.^2));
    Iyy=sum(TM.*(posx.^2+posz.^2));
    Izz=sum(TM.*(posx.^2+posy.^2));
    Ixy=-sum(TM.*posx.*posy);
    Ixz=-sum(TM.*posx.*posz);
    Iyz=-sum(TM.*posy.*posz);
    if ~isempty(TCheck.sim_file) %check if sim file is present        
        TS=squeeze(out.sig_amp_mat(:,:,:,k));
        Ixx_sim=sum(TS(:).*(YY(:).^2+ZZ(:).^2));
        Iyy_sim=sum(TS(:).*(XX(:).^2+ZZ(:).^2));
        Izz_sim=sum(TS(:).*(XX(:).^2+YY(:).^2));
        Ixy_sim=-sum(TS(:).*XX(:).*YY(:));
        Ixz_sim=-sum(TS(:).*XX(:).*ZZ(:));
        Iyz_sim=-sum(TS(:).*YY(:).*ZZ(:));
    end
    
    % combine to matrix
    TCheck.MomentsOfInertia.IMat(:,:,k)=[Ixx, Ixy, Ixz ; Ixy, Iyy, Iyz; Ixz, Iyz, Izz]/sum(TM(:));
    if ~isempty(TCheck.sim_file) %check if sim file is present    
        TCheck.MomentsOfInertia.IMat_sim(:,:,k)=[Ixx_sim, Ixy_sim, Ixz_sim ; Ixy_sim, Iyy_sim, Iyz_sim; Ixz_sim, Iyz_sim, Izz_sim]/sum(TS(:));
    end
 
    % eigenvectors and values
    [TCheck.MomentsOfInertia.V(:,:,k),TCheck.MomentsOfInertia.D(:,:,k)] = eig(TCheck.MomentsOfInertia.IMat(:,:,k));
    if ~isempty(TCheck.sim_file) %check if sim file is present  
        [TCheck.MomentsOfInertia.V_sim(:,:,k),TCheck.MomentsOfInertia.D_sim(:,:,k)] = eig(TCheck.MomentsOfInertia.IMat_sim(:,:,k));
    end
    
    % determine center of mass above threshhold
    AT=TM.*(TM/max(TM)>TCheck.CenterOfMass.thresh); 
    TCheck.CenterOfMass.cMassMat(2,k)=sum(AT.*posx)/sum(AT); %% INDEX FLIPPED
    TCheck.CenterOfMass.cMassMat(1,k)=sum(AT.*posy)/sum(AT); %% INDEX FLIPPED
    TCheck.CenterOfMass.cMassMat(3,k)=sum(AT.*posz)/sum(AT);
    
    % mean and max value above threshhold
    TCheck.CenterOfMass.cMassMean(k)=mean(AT((TM/max(TM)>TCheck.CenterOfMass.thresh)));
    TCheck.CenterOfMass.cMassMax(k)=max(AT);
    
    % determine center of mass above threshhold
    if ~isempty(TCheck.sim_file) %check if sim file is present
        AS=TS(:).*(TS(:)/max(TS(:))>TCheck.CenterOfMass.thresh);
        TCheck.CenterOfMass.cMassMat_sim(1,k)=sum(AS.*XX(:))/sum(AS);
        TCheck.CenterOfMass.cMassMat_sim(2,k)=sum(AS.*YY(:))/sum(AS);
        TCheck.CenterOfMass.cMassMat_sim(3,k)=sum(AS.*ZZ(:))/sum(AS);
    end
    
end

% determine radius from center of mass above threshhold
TCheck.CenterOfMass.RfromCMass=sqrt((TCheck.sensorPOS(:,1)*1e3-TCheck.CenterOfMass.cMassMat(1,:)').^2+(TCheck.sensorPOS(:,2)*1e3-TCheck.CenterOfMass.cMassMat(2,:)').^2);
if ~isempty(TCheck.sim_file) %check if sim file is present
    TCheck.CenterOfMass.RfromCMass_sim=sqrt((TCheck.probe.Sensors(:,1)*1e3-TCheck.CenterOfMass.cMassMat_sim(1,:)').^2+(TCheck.probe.Sensors(:,2)*1e3-TCheck.CenterOfMass.cMassMat_sim(2,:)').^2);
end

% plot
figure;
plot(TCheck.CenterOfMass.cMassMat(1,:),TCheck.CenterOfMass.cMassMat(2,:),'.'); hold;
if ~isempty(TCheck.sim_file) %check if sim file is present
    plot(TCheck.CenterOfMass.cMassMat_sim(1,:),TCheck.CenterOfMass.cMassMat_sim(2,:),'.');
end
plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'k.'); 
plot(0,0,'o'); hold;
axis image
xlim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
ylim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
if ~isempty(TCheck.sim_file) %check if sim file is present
    legend('Measurement','Simulation')
end
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off');
saveas(gcf,[TCheck.path,'\figures\orientation\center_of_mass.emf'])
close;

figure; 
pvec = TCheck.CenterOfMass.RfromCMass;
plot(pvec, '.-'); hold;
if ~isempty(TCheck.sim_file) %check if sim file is present
    plot(TCheck.CenterOfMass.RfromCMass_sim, '.-');
end
hold;
xlabel('Channel Number')
ylabel('In-Plane Detector - Focus Distance (mm)')
xlim([1, TCheck.probe.NumDetectors])
ylim([min(pvec)-0.1*max(pvec), max(pvec)+0.1*max(pvec)])
if ~isempty(TCheck.sim_file) %check if sim file is present
    legend('Measurement','Simulation')
end
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\orientation\radius.emf'])
close;
clear pvec

figure; 
pvec = TCheck.CenterOfMass.cMassMat(3,:)-mean(TCheck.CenterOfMass.cMassMat(3,:));
plot(pvec, '.-');
if ~isempty(TCheck.sim_file) %check if sim file is present
    plot(TCheck.CenterOfMass.cMassMat_sim(3,:), '.-'); 
end
xlabel('Channel Number')
ylabel('Out-of-Plane Focus Offset (mm)')
xlim([1, TCheck.probe.NumDetectors])
ylim([min(pvec)-0.1*max(pvec), max(pvec)+0.1*max(pvec)])
if ~isempty(TCheck.sim_file) %check if sim file is present
    legend('Measurement','Simulation')
end
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\orientation\out_of_plane_offset.emf'])
close;
clear pvec

figure; 
pvec1 = 180*acos(squeeze(TCheck.MomentsOfInertia.V(1,3,:)))/pi-90;
pvec2 = 180*acos(squeeze(TCheck.MomentsOfInertia.V(2,3,:)))/pi-90;
plot(pvec1,'.-'); hold;
plot(pvec2,'.-'); hold;
xlabel('Channel Number')
ylabel('Tilt Angle (°)')
xlim([1, TCheck.probe.NumDetectors])
ylim([min(min(pvec1),min(pvec2))-0.1*max(max(pvec1),max(pvec2)), max(max(pvec1),max(pvec2))+0.1*max(max(pvec1),max(pvec2))])
legend('X','Y')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\orientation\tilt_angle.emf'])
close;
% 
% figure;
% plot(squeeze(TCheck.MomentsOfInertia.IMat(1,1,:)), '.-'); hold;
% plot(squeeze(TCheck.MomentsOfInertia.IMat(2,2,:)), '.-');
% plot(squeeze(TCheck.MomentsOfInertia.IMat(3,3,:)), '.-'); 
% if ~isempty(TCheck.sim_file) %check if sim file is present
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(1,1,:)), '.-'); 
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(2,2,:)), '.-');
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(3,3,:)), '.-'); 
% end
% hold;
% xlabel('Channel Number')
% xlim([1, TCheck.probe.NumDetectors])
% legend('XX', 'YY', 'ZZ','XX_{sim}', 'YY_{sim}', 'ZZ_{sim}')
% set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
% saveas(gcf,[TCheck.path,'\figures\orientation\Imom_diagonal.emf'])
% close;
% 
% figure;
% plot(squeeze(TCheck.MomentsOfInertia.IMat(1,2,:)), '.-'); hold;
% if ~isempty(TCheck.sim_file) %check if sim file is present
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(1,2,:)), '.-');
% end
% hold;
% xlabel('Channel Number')
% xlim([1, TCheck.probe.NumDetectors])
% legend('XY', 'XY_{sim}')
% set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
% saveas(gcf,[TCheck.path,'\figures\orientation\Imom_XY.emf'])
% close;
% 
% figure;
% plot(squeeze(TCheck.MomentsOfInertia.IMat(1,3,:)), '.-'); hold;
% plot(squeeze(TCheck.MomentsOfInertia.IMat(2,3,:)), '.-'); 
% if ~isempty(TCheck.sim_file) %check if sim file is present
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(1,3,:)), '.-'); 
%     plot(squeeze(TCheck.MomentsOfInertia.IMat_sim(2,3,:)), '.-'); 
% end
% hold
% xlabel('Channel Number')
% xlim([1, TCheck.probe.NumDetectors])
% legend('XZ', 'YZ', 'XZ_{sim}', 'YZ_{sim}')
% set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
% saveas(gcf,[TCheck.path,'\figures\orientation\Imom_XZ_YZ.emf'])
% close;

% combined tilt
TCheck.MomentsOfInertia.TiltComb=sqrt((180*acos(squeeze(TCheck.MomentsOfInertia.V(1,3,:)))/pi-90).^2+(180*acos(squeeze(TCheck.MomentsOfInertia.V(2,3,:)))/pi-90).^2);


%% sensitivity maps

mkdir([TCheck.path,'\figures\smaps\'])

for j=1:TCheck.probe.NumDetectors
    
    j
    % grid vector for xy slice
    xv_slice1=linspace(-1,+1,100)*TCheck.sensorPOS(j,1)*1e3;
    yv_slice1=linspace(-1,+1,100)*TCheck.sensorPOS(j,2)*1e3; 
    
    % grid vector for vertical slice
    zv_sclice1=linspace(info.zMin_scan,info.zMax_scan,info.zStep_scan+1)-info.zCenter_scan;
    xv_slice2=linspace(-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3,100);
    yv_slice2=linspace(-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3,100);

    % meshgrids
    [Xsl1,Zsl1]=meshgrid(xv_slice1,zv_sclice1);
    [Ysl1,Zsl1]=meshgrid(yv_slice1,zv_sclice1);
    [Xsl2,Ysl2]=meshgrid(xv_slice2,yv_slice2);
    
    %xy vector to cMass
    cMassVec=TCheck.CenterOfMass.cMassMat(:,j)'-TCheck.sensorPOS(j,:)*1e3;
       
    %generate slices
    figure
    SL1=slice(X,Y,Z,squeeze(TCheck.AmpMat(:,:,:,j)),Xsl1,Ysl1,Zsl1); hold;
    SZSlice=SL1.CData;
    SL2=slice(X,Y,Z,squeeze(TCheck.AmpMat(:,:,:,j)),Xsl2,Ysl2,zeros(size(Xsl2)));
    XYSlice=SL2.CData;
    set(findobj(gca,'Type','Surface'),'EdgeColor','none')
    plot3(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,TCheck.sensorPOS(:,3)*1e3,'.')
    %quiver3(TCheck.CenterOfMass.cMassMat(1,j),TCheck.CenterOfMass.cMassMat(2,j),TCheck.CenterOfMass.cMassMat(3,j),TCheck.MomentsOfInertia.V(1,1,j)*5,TCheck.MomentsOfInertia.V(2,1,j)*5,TCheck.MomentsOfInertia.V(3,1,j)*5,'r'); 
    %quiver3(TCheck.CenterOfMass.cMassMat(1,j),TCheck.CenterOfMass.cMassMat(2,j),TCheck.CenterOfMass.cMassMat(3,j),TCheck.MomentsOfInertia.V(1,2,j)*5,TCheck.MomentsOfInertia.V(2,2,j)*5,TCheck.MomentsOfInertia.V(3,2,j)*5,'r'); 
    quiver3(TCheck.CenterOfMass.cMassMat(1,j),TCheck.CenterOfMass.cMassMat(2,j),TCheck.CenterOfMass.cMassMat(3,j),TCheck.MomentsOfInertia.V(1,3,j)*5,TCheck.MomentsOfInertia.V(2,3,j)*5,TCheck.MomentsOfInertia.V(3,3,j)*5,'r'); hold;
    axis image
    close 
    
    %color limits
    min_lim=0;
    max_lim=max(max(max(TCheck.AmpMat(:,:,:,j))));  
      
    % show xy slice   
    figure
    axes('Units','pixels','Position',[100 10 400 400]);
    imagesc(xv_slice2,yv_slice2,XYSlice, [min_lim, max_lim]); hold;
    plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'w.'); 
    plot(TCheck.sensorPOS(j,1)*1e3,TCheck.sensorPOS(j,2)*1e3,'ro'); 
    plot(TCheck.CenterOfMass.cMassMat(1,j),TCheck.CenterOfMass.cMassMat(2,j),'kx'); 
    plot(0,0,'ko');
    plot([TCheck.sensorPOS(j,1)*1e3,-TCheck.sensorPOS(j,1)*1e3],[TCheck.sensorPOS(j,2)*1e3,-TCheck.sensorPOS(j,2)*1e3],'w--')
    colormap jet
    axis image
    colorbar
    xlabel('x (mm)')
    ylabel('y (mm)')
    yylim=get(gca,'ylim');
    xxlim=get(gca,'xlim');
    myt=text(0.3*xxlim(2),0.7*yylim(2),['Focus:  ', num2str(TCheck.CenterOfMass.RfromCMass(j),'%.2f'), ' mm']);
    myt.Color='w';    
    myt=text(0.3*xxlim(2),0.85*yylim(2),['Sens_{max}:  ', num2str(TCheck.CenterOfMass.cMassMax(j),'%.0f') ]);
    myt.Color='w';
    grid on
    set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
    saveas(gcf,[TCheck.path,'\figures\smaps\XY',num2str(j),'.emf'])
    close;
    
    %show vertical slice
    figure;
    axes('Units','pixels','Position',[100 10 400 400]);    
    imagesc([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3],zv_sclice1,flip(SZSlice,2),[min_lim, max_lim]); hold
    plot(-TCheck.probe.Radius*1e3,0,'ro'); 
    plot([-TCheck.probe.Radius*1e3,TCheck.probe.Radius*1e3],[0,0],'w--')
    plot(0,0,'ko'); 
    plot(-TCheck.probe.Radius*1e3+sqrt(cMassVec(1).^2+cMassVec(2).^2),TCheck.CenterOfMass.cMassMat(3,j),'kx');
    hold;
    colormap jet
    cb=colorbar;
    cb.Visible='off';
    axis image
    xlabel('Section Position (mm)')
    ylabel('z (mm)')
    yylim=get(gca,'ylim');
    xxlim=get(gca,'xlim');
    myt=text(0.3*xxlim(2),0.6*yylim(2),['Z-Offset:  ', num2str(TCheck.CenterOfMass.cMassMat(3,j),'%.2f'), ' mm' ]);
    myt.Color='w';
    set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
    saveas(gcf,[TCheck.path,'\figures\smaps\XZ',num2str(j),'.emf'])
    close;
    
    
end

% show xy mean slice   

%color limits
min_lim=0;
max_lim=max(max(max(mean(TCheck.AmpMat(:,:,maxSliceNum,:),4))));  

figure
imagesc([min(posx),max(posx)],[min(posy),max(posy)],interp2(mean(TCheck.AmpMat(:,:,maxSliceNum,:),4),2), [min_lim, max_lim]); hold;
%imagesc([min(posx),max(posx)],[min(posy),max(posy)],interp2(flip(mean(TCheck.AmpMat(:,:,maxSliceNum,:),4),1),2), [min_lim, max_lim]); hold;
plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'k.'); 
plot(0,0,'ko'); hold;
colormap jet
axis image
colorbar
xlabel('x (mm)')
ylabel('y (mm)')
xlim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
ylim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\smaps\XYmean.emf'])
close;

%% elevational resolution

% generate output path
mkdir([TCheck.path,'\figures\resolution\'])

% generate resolution maps
TCheck.Resolution.ElResMat=zeros(size(TCheck.AmpMat,1),size(TCheck.AmpMat,2),TCheck.probe.NumDetectors);
TCheck.Resolution.ElResMatSim=zeros(size(out.sig_amp_mat,1),size(out.sig_amp_mat,2),TCheck.probe.NumDetectors);
TCheck.Resolution.ElResTot=zeros(size(TCheck.AmpMat,1),size(TCheck.AmpMat,2));
TCheck.Resolution.ElResTotSim=zeros(size(out.sig_amp_mat,1),size(out.sig_amp_mat,2));

for k=1:TCheck.probe.NumDetectors
    k
    
    % measured resolution
    for i=1:size(TCheck.AmpMat,1)
        for j=1:size(TCheck.AmpMat,2)  
           try TCheck.Resolution.ElResMat(i,j,k)=get_fwhm(linspace(info.zMin_scan,info.zMax_scan,info.zStep_scan)-info.zCenter_scan,squeeze(TCheck.AmpMat(i,j,:,k))); end  
        end
    end
    if ~isempty(TCheck.sim_file) %check if sim file is present
        % simulated
        for i=1:size(out.sig_amp_mat,1)
            for j=1:size(out.sig_amp_mat,2)
               try TCheck.Resolution.ElResMatSim(i,j,k)=get_fwhm(poszz,squeeze(out.sig_amp_mat(i,j,:,k)));  end
            end
        end
    end
end

% total measured
for i=1:size(TCheck.AmpMat,1)
    for j=1:size(TCheck.AmpMat,2)
        TCheck.Resolution.ElResTot(i,j)=get_fwhm(linspace(info.zMin_scan,info.zMax_scan,info.zStep_scan)-info.zCenter_scan,squeeze(mean(TCheck.AmpMat(i,j,:,:),4)));
    end
end

% total simulated
if ~isempty(TCheck.sim_file) %check if sim file is present
    for i=1:size(out.sig_amp_mat,1)
        for j=1:size(out.sig_amp_mat,2)
            TCheck.Resolution.ElResTotSim(i,j)=get_fwhm(poszz,squeeze(mean(out.sig_amp_mat(i,j,:,:),4)));
        end
    end
end

% plot
for j=1:TCheck.probe.NumDetectors
    tic
    j
    %h = figure('Visible','off')
    figure
    imagesc([min(posx),max(posx)],[min(posy),max(posy)],interp2(TCheck.Resolution.ElResMat(:,:,j),2),[0.5,3]); hold;
    plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'k.'); 
    plot(0,0,'ko');
    plot([TCheck.sensorPOS(j,1)*1e3,-TCheck.sensorPOS(j,1)*1e3],[TCheck.sensorPOS(j,2)*1e3,-TCheck.sensorPOS(j,2)*1e3],'k--'); hold; 

    xlim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]);
    ylim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]);
    colormap jet
    c = colorbar;
    c.Label.String = 'd (mm)';
    xlabel('x (mm)')
    ylabel('y (mm)')
    axis square
    set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
    saveas(gcf,[TCheck.path,'\figures\resolution\elevation_res',num2str(j),'.emf'])
    close;
    toc
    
end
   
figure
imagesc([min(posx),max(posx)],[min(posy),max(posy)],interp2(TCheck.Resolution.ElResTot,2),[0.5,3]); hold; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'k.'); 
plot(0,0,'ko'); hold;
colormap jet
axis image
c = colorbar;
c.Label.String = 'El. Resolution (mm)';
xlabel('x (mm)')
ylabel('y (mm)')
xlim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
ylim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\resolution\elevation_res_tot.emf'])
close;

if ~isempty(TCheck.sim_file) %check if sim file is present
    figure
    imagesc([min(posx),max(posx)],[min(posy),max(posy)],TCheck.Resolution.ElResTotSim,[1,3]); hold;
    plot(TCheck.sensorPOS(:,1)*1e3,TCheck.sensorPOS(:,2)*1e3,'k.'); 
    plot(0,0,'ko'); hold;
    colormap jet
    axis image
    colorbar
    xlabel('x (mm)')
    ylabel('y (mm)')
    xlim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
    ylim([-TCheck.probe.Radius*1e3,+TCheck.probe.Radius*1e3]*1.05);
    grid on
    set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
    saveas(gcf,[TCheck.path,'\figures\resolution\elevation_res_tot_sim.emf'])
    close;
end

for j=1:TCheck.probe.NumDetectors
    j
    TCheck.CenterOfMass.cMassMat(:,j)
    [Xcms,Ycms,Zcms]=meshgrid(...
        TCheck.CenterOfMass.cMassMat(2,j),...
        TCheck.CenterOfMass.cMassMat(1,j),...
        linspace(info.zMin_scan,info.zMax_scan,200)-info.zCenter_scan...
        );
    F = scatteredInterpolant(posx',posy',posz',double(TCheck.AmpVec(j,:)'));
    try
        TCheck.Resolution.ElRescMass(j)=get_fwhm(linspace(info.zMin_scan,info.zMax_scan,200)-info.zCenter_scan,squeeze(F(Xcms,Ycms,Zcms))); 
    catch
        'FWHM calculation failed'
    end
end

if ~isempty(TCheck.sim_file) %check if sim file is present
    for j=1:TCheck.probe.NumDetectors
        j
        TCheck.CenterOfMass.cMassMat_sim(:,j)
        [Xcms,Ycms,Zcms]=meshgrid(...
            TCheck.CenterOfMass.cMassMat_sim(1,j),...
            TCheck.CenterOfMass.cMassMat_sim(2,j),...
            linspace(info.zMin_scan,info.zMax_scan,200)-info.zCenter_scan...
            );
        F0=out.sig_amp_mat(:,:,:,j);
        F = scatteredInterpolant(XX(:),YY(:),ZZ(:),F0(:));
        F.Method = 'natural';
        TCheck.Resolution.ElRescMass_sim(j)=get_fwhm(linspace(info.zMin_scan,info.zMax_scan,200)-info.zCenter_scan,squeeze(F(Xcms,Ycms,Zcms)));
    end
end

figure;
pvec1 = squeeze(TCheck.Resolution.ElResMat(round(size(TCheck.AmpMat,1)/2),round(size(TCheck.AmpMat,2)/2),:));
pvec2 = TCheck.Resolution.ElRescMass;
plot(pvec1, '.-'); hold;
plot(pvec2, '.-');
if ~isempty(TCheck.sim_file) %check if sim file is present
    plot(squeeze(TCheck.Resolution.ElResMatSim(round(size(out.sig_amp_mat,1)/2),round(size(out.sig_amp_mat,2)/2),:)), '.-');
    plot(TCheck.Resolution.ElRescMass_sim, '.-');
end
hold;
xlim([1 TCheck.probe.NumDetectors])
ylim([min(min(pvec1),min(pvec2))-0.1*max(max(pvec1),max(pvec2)), max(max(pvec1),max(pvec2))+0.1*max(max(pvec1),max(pvec2))])
xlabel('Channel Number')
ylabel('Elevation Resolution (mm)')
legend('COR','CMS','COR_{sim}','CMS_{sim}')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\resolution\elevation_res_COR.emf'])
close;

%% check FFT spectra

% set up fft
f=[1:2030]'*(40/2030);
tvec=[1:2030]/40;

% windowing and plotting parameters
TCheck.FFTSpectra.winsize=50;
TCheck.FFTSpectra.plotSamplesHalf=200;

% generate output path
mkdir([TCheck.path,'\figures\fft\'])

TCheck.FFTSpectra.sigMatFilt_XYZcenter_Win=zeros(size(TCheck.sigMatFilt_XYZcenter));
TCheck.FFTSpectra.FsigMatFilt_XYZcenter=zeros(size(TCheck.sigMatFilt_XYZcenter));

for i=1:TCheck.probe.NumDetectors
       
    % select signal
    sig=TCheck.sigMatFilt_XYZcenter(:,i);
        
    % find peak position
    peak_pos_max=find(sig==max(sig));
    peak_pos_min=find(sig==min(sig));
    TCheck.FFTSpectra.peak_pos(i)=int16(round((peak_pos_max+peak_pos_min)/2));
    
    % apply gaussian window
    win=circshift(gausswin(2030,TCheck.FFTSpectra.winsize),TCheck.FFTSpectra.peak_pos(i)-1015);
    sig=sig.*win;
    TCheck.FFTSpectra.sigMatFilt_XYZcenter_Win(:,i)=sig;
    
    % signal maximum
    TCheck.FFTSpectra.SigMax(i)=max(sig)-min(sig);
    
    % calculate fft
    fsig=fft(sig);
    TCheck.FFTSpectra.FsigMatFilt_XYZcenter(:,i)=fsig;
    
    % determine dominant frequency using center of mass
    TCheck.FFTSpectra.fCMass(i)=sum(f(1:end/4).*abs(fsig(1:end/4)))/sum(abs(fsig(1:end/4)));
    
    % calculate bandwidth
    TCheck.FFTSpectra.BW_total(i)=get_fwhm(f,abs(fsig));
    TCheck.FFTSpectra.BW_perc(i)=TCheck.FFTSpectra.BW_total(i)./TCheck.FFTSpectra.fCMass(i);
    
    %  maximum of spectrum
    TCheck.FFTSpectra.SpecMax(i)=max(abs(fsig(1:end/4)));
    cmaxT=f(find(TCheck.FFTSpectra.SpecMax(i)==abs(fsig(1:end/4))));
    TCheck.FFTSpectra.fCMax(i)=cmaxT(1);

    %plot 
    figure;
    subplot(2,1,1)
    plot(tvec(TCheck.FFTSpectra.peak_pos(i)-TCheck.FFTSpectra.plotSamplesHalf:TCheck.FFTSpectra.peak_pos(i)+TCheck.FFTSpectra.plotSamplesHalf),sig(TCheck.FFTSpectra.peak_pos(i)-TCheck.FFTSpectra.plotSamplesHalf:TCheck.FFTSpectra.peak_pos(i)+TCheck.FFTSpectra.plotSamplesHalf)); 
    %hold; plot(tvec(TCheck.FFTSpectra.peak_pos(i)-TCheck.FFTSpectra.plotSamplesHalf:TCheck.FFTSpectra.peak_pos(i)+TCheck.FFTSpectra.plotSamplesHalf),max(sig(:)).*win(TCheck.FFTSpectra.peak_pos(i)-TCheck.FFTSpectra.plotSamplesHalf:TCheck.FFTSpectra.peak_pos(i)+TCheck.FFTSpectra.plotSamplesHalf)); hold;
    axis tight
    xlabel('Time (µs)')
    ylabel('Signal (a.u.)')
    subplot(2,1,2)
    plot(f(1:end/4), abs(fsig(1:end/4)))
    xlabel('Frequency (MHz)')
    ylabel('Magnitude (a.u.)')  
    yylim=get(gca,'ylim');
    xxlim=get(gca,'xlim');
    text(0.75*xxlim(2),0.75*yylim(2),['f_{cms}: ', num2str(TCheck.FFTSpectra.fCMass(i),'%.2f'), ' MHz'])
    text(0.75*xxlim(2),0.55*yylim(2),['BW_{%}: ', num2str(TCheck.FFTSpectra.BW_perc(i)*100,'%.2f') ])    
    set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
    set(gcf, 'Visible', 'on')
    saveas(gcf,[TCheck.path,'\figures\fft\fft_det',num2str(i),'.emf'])
    close;
    
end

figure;
pvec1 = TCheck.FFTSpectra.fCMass;
pvec2 = TCheck.FFTSpectra.fCMax;
plot(pvec1,'.-'); hold;
plot(pvec2,'.');
plot([1 TCheck.probe.NumDetectors],[TCheck.probe.CenterFrequency,TCheck.probe.CenterFrequency]*1e-6,'k'); hold;
ylim([TCheck.probe.CenterFrequency*1e-6*0.75,TCheck.probe.CenterFrequency*1e-6*1.25])
xlim([1 TCheck.probe.NumDetectors])
xlabel('Channel Number')
ylabel('Frequency (MHz)')
legend('f_{cms}','f_{max}', 'target')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\fft\f_center.emf'])
close;

figure;
plot(TCheck.FFTSpectra.SpecMax,'.-');
xlim([1 TCheck.probe.NumDetectors])
xlabel('Channel Number')
ylabel('Maximum of Spectrum (a.u.)')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\fft\spec_max.emf'])
close;

figure;
pvec1 = TCheck.FFTSpectra.SigMax;
pvec2 = TCheck.CenterOfMass.cMassMax;
plot(pvec1,'.-'); hold;
plot(pvec2, '.-'); hold;
xlim([1 TCheck.probe.NumDetectors])
ylim([min(min(pvec1),min(pvec2))-0.1*max(max(pvec1),max(pvec2)), max(max(pvec1),max(pvec2))+0.1*max(max(pvec1),max(pvec2))])
legend('COR','CMS')
xlabel('Channel Number')
ylabel('Signal Amplitude(a.u.)')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\fft\sig_max.emf'])
close;

figure;
pvec1 = TCheck.FFTSpectra.BW_perc*100;
pvec2 = ones(TCheck.probe.NumDetectors,1)*sqrt(2)*TCheck.probe.Bandwidth*100;
plot(pvec1,'.-'); hold;
plot(pvec2 ,'k')
xlim([1 TCheck.probe.NumDetectors])
%ylim([min(min(pvec1),min(pvec2))-0.1*max(max(pvec1),max(pvec2)), max(max(pvec1),max(pvec2))+0.1*max(max(pvec1),max(pvec2))])
ylim([sqrt(2)*TCheck.probe.Bandwidth*100*0.8 sqrt(2)*TCheck.probe.Bandwidth*100*1.5])
xlabel('Channel Number')
ylabel('Bandwidth (%)')
legend('measured','target')
grid on
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\fft\bandwidth.emf'])
close;


%% image reconstruction

% average center Z plane signals
ReconSigMat=mean(TCheck.sigMatFilt_Zcenter,3);

% speed of sound
c=water_sos(TCheck.Temp)+15;

% fov limits
limits(1) = TCheck.probe.Radius - TCheck.Recon.fov_size/2; 
limits(2) = TCheck.probe.Radius + TCheck.Recon.fov_size/2; 

% filter limits
ampmin=0+(TCheck.Recon.fmin/20)*1i; 
ampmax=0+(TCheck.Recon.fmax/20)*1i; 
freqVec = [0, TCheck.Recon.fmin*0.8, TCheck.Recon.fmin,     TCheck.Recon.fmax,    TCheck.Recon.fmax*1.2, 20]/20;   
ampVec =  [0, 0,         ampmin,  ampmax, 0          0]; 

% setup reconstruction
clfile = which('itheract.cl'); 
[clpath,~] = fileparts(clfile);
startupcl('vendor','amd','type','gpu','clpath',clpath);
obj = clBackprojection('2d'); 
obj.RoI = limits;
obj.N = TCheck.Recon.n;
obj.Projections = TCheck.probe.NumDetectors*2+1;
obj.Fs = 40e6;
obj.TimeResolution = 3;          
obj.Sensor = {TCheck.sensorPOS};
obj.ImpulseResponse = [];
obj.FrequencyVector = freqVec;
obj.AmplitudeVector = ampVec;
obj.Setup; 

% run recon
TCheck.Recon.R = rot90(obj.Run(uint16(ReconSigMat), 1.0, c, 0),3);

% delet and clear object
delete(obj);
clear obj;
cleanupcl;

figure;
imagesc([-1e3*(TCheck.Recon.fov_size/2), 1e3*(TCheck.Recon.fov_size/2)],[-1e3*(TCheck.Recon.fov_size/2), 1e3*(TCheck.Recon.fov_size/2)],max(TCheck.Recon.R,0))
colormap jet
axis image
grid on
xlabel('x (mm)')
ylabel('y (mm)')
set(gcf, 'Color', 'none', 'Inverthardcopy', 'off'); 
saveas(gcf,[TCheck.path,'\figures\recon_center_plane.emf'])
%close;

%% dead elements (peak -6dB below mean)

TCheck.DeadElements=sum(TCheck.CenterOfMass.cMassMax<10^(-6/10)*mean(TCheck.CenterOfMass.cMassMax));

%% measured radius

TCheck.RadiusMeasMean=mean(1e3*water_sos(TCheck.Temp)*double(TCheck.FFTSpectra.peak_pos)/4e7);
TCheck.RadiusMeasStd=std(1e3*water_sos(TCheck.Temp)*double(TCheck.FFTSpectra.peak_pos)/4e7);

%% create power point report

% Open existing presentation
isOpen  = exportToPPTX();
if ~isempty(isOpen),
    % If PowerPoint already started, then close first and then open a new one
    exportToPPTX('close');
end
exportToPPTX('open','templateB.pptx');

% add title slide
exportToPPTX('addslide','Master',1,'Layout','Titelfolie');
exportToPPTX('addtext','Transducer Test Report','Position','Titel');
exportToPPTX('addtext',['Design ID: ', TCheck.UTDID, '    SN: ', TCheck.UTSN ],'Position','Untertitel');
exportToPPTX('addtext',datestr(now,'mmmm dd, yyyy HH:MM:SS'),'Position','Text');

% add summary slide 
exportToPPTX('addslide','Master',1,'Layout','Summary');
exportToPPTX('addtext','Summary','Position','Titel');

exportToPPTX('addtext',...
    sprintf([...
    'Center frequency:    \r',...
    'Bandwidth (FWHM):     \r',...
    'Sensitivity variation:    \r',...
    'Functional elements:    \r',...
    'Radius:   \r',...
    'Focus:    \r',...
    'Tilt:    \r',...
    'Elev. offset:    \r',...
    'Elev. resolution: \r',...
    '\r',...
    'Angular coverage:    \r',...
    'Curvature:   \r',...
    'El. height:    \r',...,
    'Pitch:    \r'...
    ]),...
    'Position','Inhaltsplatzhalter 1');


exportToPPTX('addtext',...
    sprintf([...
...
    num2str(mean(TCheck.FFTSpectra.fCMass),'%.1f'), ' +- ', num2str(std(TCheck.FFTSpectra.fCMass),'%.1f') , ' MHz (target: ', num2str(TCheck.probe.CenterFrequency/1e6) ,' MHz) \r',... % center frequency
...
    num2str(mean(TCheck.FFTSpectra.BW_perc*100),'%.0f'), ' +- ', num2str(std(TCheck.FFTSpectra.BW_perc*100),'%.0f') ,' %% (target: ', num2str(sqrt(2)*TCheck.probe.Bandwidth*100,'%.0f'),' %%)  \r',... % bandwidth
... 
    num2str(20*log10(max(TCheck.CenterOfMass.cMassMax/mean(TCheck.CenterOfMass.cMassMax))),'%.0f')  ,' dB \r',... % sensitivty variation in cms
...   
    num2str(TCheck.probe.NumDetectors-TCheck.DeadElements), ' out of ' num2str(TCheck.probe.NumDetectors),'   \r',... %functional elements
...
    num2str(TCheck.RadiusMeasMean,'%.1f'),' +- ',num2str(TCheck.RadiusMeasStd,'%.1f'),' mm (target: ',num2str(TCheck.probe.Radius*1e3),' mm) \r',... %radius
...    
    num2str(mean(TCheck.CenterOfMass.RfromCMass), '%.1f'),' +- ', num2str(std(TCheck.CenterOfMass.RfromCMass), '%.1f'),' mm \r',... %average focus
...
    num2str(mean(TCheck.MomentsOfInertia.TiltComb), '%.1f'),' +- ', num2str(std(TCheck.MomentsOfInertia.TiltComb), '%.1f'),' ° \r',... % average tilt
...
    'max: ',num2str(max(TCheck.CenterOfMass.cMassMat(3,:)-mean(TCheck.CenterOfMass.cMassMat(3,:))), '%.2f'),' mm / min:', num2str(min(TCheck.CenterOfMass.cMassMat(3,:)-mean(TCheck.CenterOfMass.cMassMat(3,:))), '%.2f'),' mm \r',... % average elev. offset
...
    num2str(mean(TCheck.Resolution.ElRescMass), '%.2f'),' +- ', num2str(std(TCheck.Resolution.ElRescMass), '%.2f'),' mm (target: ','-', ' mm) \r',... %elevational resolution
...
    '\r',...
...
    num2str(TCheck.probe.Coverage),'°  \r',... %angular coverage
...    
    num2str(TCheck.probe.Curvature*1e3),' mm  \r',... %curvature
...
    num2str(round(TCheck.probe.Height*1e3,2)),' mm  \r',... %height
...
    num2str(round(TCheck.probe.Pitch*1e3,3)),' mm  \r',... %pitch
...
    ]),...
    'Position','Inhaltsplatzhalter 2');

% add frequency slide 
exportToPPTX('addslide','Master',1,'Layout','Frequency');
exportToPPTX('addtext','Frequency and Bandwidth','Position','Titel');
exportToPPTX('addtext','Center frequency','Position','Text Placeholder 1');
exportToPPTX('addtext','Bandwidth','Position','Text Placeholder 2');
exportToPPTX('addtext','Amplitude (peak-to-peak)','Position','Text Placeholder 3');
exportToPPTX('addtext',sprintf(['Mean (CMS):' ,'\r',  '', num2str(mean(TCheck.FFTSpectra.fCMass),'%.1f'), ' +- ', num2str(std(TCheck.FFTSpectra.fCMass),'%.1f') ,' MHz']),'Position','Text Placeholder 4');
exportToPPTX('addtext',sprintf(['Mean:' ,'\r',  '', num2str(mean(TCheck.FFTSpectra.BW_perc*100),'%.0f'), ' +- ', num2str(std(TCheck.FFTSpectra.BW_perc*100),'%.0f') ,' %%']),'Position','Text Placeholder 5');
exportToPPTX('addtext',sprintf(['Mean (Total):' ,'\r',  '', num2str(mean(TCheck.CenterOfMass.cMassMax),'%.0f'), ' +- ', num2str(std(TCheck.CenterOfMass.cMassMax),'%.0f') ,' ']),'Position','Text Placeholder 6');
exportToPPTX('addpicture',[TCheck.path,'\figures\fft\f_center.emf'],'Position','centerf','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\fft\bandwidth.emf'],'Position','bandwidth','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\fft\sig_max.emf'],'Position','amplitude','Scale','max');

% add sensitivity slide
exportToPPTX('addslide','Master',1,'Layout','Sensitivity');
exportToPPTX('addtext','Sensitivity and Resolution','Position','Titel');
exportToPPTX('addtext','Total Sensitivity','Position','Text Placeholder 1');
exportToPPTX('addtext','Out-of-plane Resolution','Position','Text Placeholder 2');
exportToPPTX('addpicture',[TCheck.path,'\figures\smaps\XYmean.emf'],'Position','senstot','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\resolution\elevation_res_tot.emf'],'Position','elres','Scale','max');

% add orientation slide
exportToPPTX('addslide','Master',1,'Layout','Orientation');
exportToPPTX('addtext','Orientation and Focus','Position','Titel');
exportToPPTX('addpicture',[TCheck.path,'\figures\orientation\center_of_mass.emf'],'Position','cmass','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\orientation\radius.emf'],'Position','focus','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\orientation\out_of_plane_offset.emf'],'Position','zoff','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\orientation\tilt_angle.emf'],'Position','tilt','Scale','max');
exportToPPTX('addpicture',[TCheck.path,'\figures\resolution\elevation_res_COR.emf'],'Position','elrescmass','Scale','max');

for i=1:TCheck.probe.NumDetectors
 
    i
   
    % add element slide
    exportToPPTX('addslide','Master',1,'Layout','Element');
    exportToPPTX('addtext',['Element No. ', num2str(i)],'Position','Titel');
    exportToPPTX('addtext','Signal at COR','Position','Text Placeholder 1');
    exportToPPTX('addtext','Sensitivity','Position','Text Placeholder 2');
    exportToPPTX('addtext','Slice Thickness','Position','Text Placeholder 3');
    exportToPPTX('addpicture',[TCheck.path,'\figures\fft\fft_det',num2str(i),'.emf'],'Position','signal','Scale','max');
    exportToPPTX('addpicture',[TCheck.path,'\figures\smaps\XY',num2str(i),'.emf'],'Position','sensxy','Scale','max');
    exportToPPTX('addpicture',[TCheck.path,'\figures\smaps\XZ',num2str(i),'.emf'],'Position','senssec','Scale','max');
    exportToPPTX('addpicture',[TCheck.path,'\figures\resolution\elevation_res',num2str(i),'.emf'],'Position','slicet','Scale','max');
    
end


% save
newFile = exportToPPTX('save', ['Transducer Test Report ',TCheck.UTDID,' ',TCheck.UTSN,' (',datestr(now,'yyyy-mm-dd-HH-MM-SS'),')']);

% close presentation (and clear all temporary files)
exportToPPTX('close');

fprintf('New file has been saved: <a href="matlab:winopen(''%s'')">%s</a>\n',newFile,newFile);

%% TO DO

% slide with 3 reconstructed images at z-1,0,1 mm
% slide with simulation

 %num2str(mean(TCheck.CenterOfMass.RfromCMass), '%.0f'),' +- ', num2str(std(TCheck.CenterOfMass.RfromCMass), '%.0f'),' mm (target: ', num2str(mean(TCheck.CenterOfMass.RfromCMass_sim), '%.0f'), ' mm) \r',... %average focus
