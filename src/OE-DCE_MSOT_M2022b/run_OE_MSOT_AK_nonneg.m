function [segm] = run_OE_MSOT_AK_nonneg (msot_file_medair, msot_file_02, saving_DIR)

% Add MSOT MAtlab library to path
% This library was created by iThera
%addpath( genpath('com.itheramedical.msotlib_beta_rev723') )
%javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\xbean.jar %%AK
%javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\msotbeans.jar %%AK

%% loadMSOTMsp  Load MSOT Specrally Unmixed Data
%   [R wls zpos ts datainfo] = LOADMSOTMSP()
%           opens Open File Dialog to choose file, then loads the first
%           Msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(file)
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo)
%           opens the specified file, then loads the first msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum)
%           opens the specified file, then loads the specified msp
%   [R wls zpos ts datainfo] = LOADMSOTMSP(datainfo,mspNum,selMat)
%           opens the specified file, then loads the specified msp
%           restricted to selMat, which is a multidimensional structure of
%           frame ids. For all, use datainfo.MspNode(i).Structure
%
% This function loads the MSOT MSPs from the selected scan folder.
%
% Return values:
% R          An array containing the MSPed images in 8D
%            1) x
%            2) y
%            3) z (only if 3D detector)
%            4) RUN (repetition of the whole acquisition step, time
%               dimension (timepoints as separate vector)
%            5) z-position (if 2D system with translation stage)
%            6) REPETITION (currently unused)
%            7) Wavelength (see wls vector)
%            8) Individual Images (only if not averaged)
% wls        Wavelength Vector
% zpos       Array of z-stage positions
% ts         Multidimensional array of timestamps in s, same
%            dimension as R
%%

% Constant
linewidthTK = 2;
roi_approve = 0;

for idx =1:2
    
    switch idx
        case 1
            msot_file = msot_file_medair;
        case 2
            msot_file = msot_file_02;
    end
    datainfo = loadMSOT(msot_file);
    
    % this loads the MSP files
    [R, ~, ts, ~] = loadMSOTMsp(msot_file);
    msot = squeeze(R);
    time_min = (squeeze( ts(:,:,:,1) ) / 60)';
    
    % obtain the number of spectrally unmixed images
    [d1,d2,~,~,n_tmps,n_msp] = size(R);
    
    % Loads the RECON files
    [R_recon, ~] = loadMSOTRecon(datainfo);
    msot_recon = squeeze(R_recon);
    ref_anat = msot_recon(:,:,1,5); % image before injection at 875 nm, first timepoint
    
    
    %% Define ROIs and approve
    
    % Use Oxy-hemoglobin signal for ROI identification
    deoxy = squeeze(msot(:,:,:,1));
    oxy = squeeze( msot(:,:,:,n_msp));
    total = zeros(d1,d2,n_tmps);
    
    % Assign negative values to zero
    deoxy(deoxy < 0) = 0;
    oxy(oxy < 0 ) = 0;
    
    for t=1:n_tmps
        total(:,:,t) = oxy(:,:,t) + deoxy(:,:,t);
    end
    
    while roi_approve == 0
        
        ref = oxy(:,:,1);
        
        %TOI
        %         ref2 = flip(ref,2); %flip left/right
        %         ref3 = flip(ref2,1); %flip top/bottom
        %         ref4 = imadjust(ref3); %Adjust intensity values
        %         figure(1); imshow(ref, 'Colormap', gray); title('Draw tumor'); %h=drawfreehand;
        %         figure(2); imshow(ref4, 'Colormap', gray); title('Draw tumor'); %h=drawfreehand;
        figure(3); imagesc(ref); colormap(gray); axis image off; title('Draw tumor');
        set(gca, 'ydir', 'normal');
        set(gca, 'xdir', 'reverse');
        h=drawfreehand;
        segm = createMask(h);
        
        %         figure(5); imshow(segm, 'Colormap', gray); title('Draw tumor'); %h=drawfreehand;
        
        str = input('Are you satisfied with the ROI? Y/N [Y]: ','s');
        if strcmpi(str,'Y')
            roi_approve = 1;
        end
        
        close all;
    end
    
    figure(4); imagesc(segm); colormap(gray); axis image; title('Draw tumor');
    set(gca, 'ydir', 'normal');
    set(gca, 'xdir', 'reverse');
    
    %% Directory for saving results
    cd(saving_DIR);
    filename=[datainfo.Name,'_OE_MSOT_MSP'];
    mkdir(filename);
    cd(filename);
    
    %% OE analysis
    map = zeros(d1,d2,n_tmps);
    
    for t = 1:n_tmps
        map(:,:,t) = oxy(:,:,t)./total(:,:,t);
    end
    
    % if total hemoglobin is 0 then some voxels in the map will be NaN. Use the
    % following line to replace NaN with zero.
    map(isnan(map))=0;
    
    map = map .*segm;
    
    percent = zeros(1,n_tmps);
    t_oxy = zeros(1,n_tmps);
    t_deoxy = zeros(1,n_tmps);
    t_total = zeros(1,n_tmps);
    
    for i=1:n_tmps
        for j=1:d1
            for k=1:d2
                if segm(j,k)>0
                    t_oxy(1,i) = t_oxy(1,i)+oxy(j,k,i);
                    t_deoxy(1,i) = t_deoxy(1,i)+deoxy(j,k,i);
                    t_total(1,i) = t_total(1,i)+total(j,k,i);
                end
            end
        end
        percent(1,i)=t_oxy(1,i)./t_total(1,i);
    end
    
    
    %% Plot average signal for oxy, deoxy, total and percent within ROI
    disp(['Average Oxy        : ' num2str(mean(t_oxy))])
    disp(['Average Deoxy      : ' num2str(mean(t_deoxy))])
    disp(['Average Oxy + Deoxy: ' num2str(mean(t_total))])
    
    figure (1);
    hold on;
    plot(time_min, t_oxy, 'bo-', 'LineWidth',linewidthTK);
    plot(time_min, t_deoxy, 'ro-', 'LineWidth',linewidthTK);
    plot(time_min, t_total, 'go-', 'LineWidth',linewidthTK);
    hold off;
    xlabel('Time (min)');
    ylabel('Signal ');
    legend ({'Oxy','Deoxy','Total'},'Location','northeast');
    savefig('Oxy_Deoxy_Total.fig');
    %     saveas(gcf, 'Oxy_Deoxy_Total.tiff');
    print(gcf, '-dtiff', 'Oxy_Deoxy_Total.tiff');
    
    figure(2)
    plot(time_min, percent, 'co-', 'LineWidth',linewidthTK);
    xlabel('Time (min)');
    ylabel('Oxy/Total');
    savefig('Oxy_over_Total.fig');
    %     saveas(gcf, 'Oxy_over_Total.tiff');
    print(gcf, '-dtiff', 'Oxy_over_Total.tiff');
    
    
    average_map = mean(map,3);
    
    figure(3)
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, average_map, 'alphadata', segm>0);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('S02');
    colorbar;
    savefig('Average S02.fig');
    %     saveas(gcf, 'Average S02.tiff');
    print(gcf, '-dtiff', 'Average S02.tiff');
    
    %     disp("Average Oxy:")
    average_oxy = mean(oxy,3);
    
    
    figure(4)
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, average_oxy, 'alphadata', segm>0);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('S02');
    colorbar;
    savefig('Average oxy.fig');
    %     saveas(gcf, 'Average oxy.tiff');
    print(gcf, '-dtiff', 'Average oxy.tiff');
    
    %     disp("Average Deoxy:")
    average_deoxy = mean(deoxy,3);
    
    figure(5)
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, average_deoxy, 'alphadata', segm>0);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('S02');
    colorbar;
    savefig('Average deoxy.fig');
    %     saveas(gcf, 'Average deoxy.tiff');
    print(gcf, '-dtiff', 'Average deoxy.tiff');
    
    
    close all;
    
    namefile = ['Result_',filename,'.mat'];
    
    switch idx
        case 1
            map_medair = average_map;
        case 2
            map_02 = average_map;
    end
    
    save(namefile, 't_oxy', 't_deoxy', 't_total','percent', 'map', 'segm','oxy','deoxy', 'ref_anat', 'time_min');
    
    
end

delta = zeros(d1,d2);
delta = map_02 - map_medair;

figure;
ax1 = axes;
imagesc(ref_anat)
axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
colormap(ax1, 'gray');
ax2 = axes;
imagesc(ax2, delta, 'alphadata', segm>0);
axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
colormap(ax2,'jet');
caxis(ax2, [-1 1]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
title('deltaS02');
colorbar;
savefig('deltaS02.fig');
% saveas(gcf, 'deltaS02.tiff');
print(gcf, '-dtiff', 'deltaS02.tiff');

delta_avg = mean(nonzeros(delta));
save('deltaS02.mat', 'map_02', 'map_medair', 'ref_anat', 'delta', 'delta_avg');

close all;


