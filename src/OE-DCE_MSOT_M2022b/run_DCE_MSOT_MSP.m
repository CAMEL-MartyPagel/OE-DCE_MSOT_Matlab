function [] = run_DCE_MSOT_MSP (msot_file, saving_DIR, segm)

% Add MSOT MAtlab library to path 
% This library was created by iThera
addpath( genpath('com.itheramedical.msotlib_beta_rev723') )
javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\xbean.jar %%AK
javaaddpath com.itheramedical.msotlib_beta_rev723\MSOTBeans\msotbeans.jar %%AK

%%
%   loadMSOTMsp  Load MSOT Specrally Unmixed Data
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

if ~exist('segm','var')
     % third parameter does not exist, so default it to something
      segm = [];
      toi_exist = 0;
else
    toi_exist = 1;
    mask_toi = segm;
end

datainfo = loadMSOT(msot_file);

% Constant
linewidthTK = 2;
fontsz = 11.5;
map_no = 4;

% this loads the entire experiment
[R, ~, ts, ~] = loadMSOTMsp(msot_file);
msot = squeeze(R);

% obtain the number of spectrally unmixed images 
[d1,d2,d3,d4,d5,n_msp] = size(R);

%% Define ROIs and approve
roi_approve = 0;

% Use ICG signal for ROI identification
prompt = 'Select chromophore for segmentation: ';
str = input(prompt,'s');
cl = str2double(str);
I = squeeze( msot(:,:,:,cl)); 

%------------------- Get Signals from all voxels ------------------------%
% [rows,cols, timepoints ] = size(I);
% dd = 1;
% dd_max = rows * cols * timepoints;
% disp(dd_max)
% 
% voxel_row = zeros(1, dd_max);
% voxel_col = zeros(1, dd_max);
% signal = zeros(dd_max,timepoints);
% 
% 
% for r = 1:rows
%     for c = 1:cols
%         for t = 1:timepoints
%             voxel_row(dd) = r;
%             voxel_col(dd) = c;
%             signal(dd, t) = I(r, c, t);
%             dd = dd + 1;
%         end
%     end
% end
% 
% disp("size(signal) = "); disp(size(signal));
% disp("size(voxel_row) = "); disp(size(voxel_row));
% disp("size(voxel_col) = "); disp(size(voxel_col));
% disp ("signal(1:10,1:100) ="); disp (signal(1:10,1:100));
% disp ("voxel_row(1,1:100) ="); disp (voxel_row(1,1:100));
% disp ("voxel_col(1,1:100) ="); disp (voxel_col(1,1:100));

%------------- Use deoxy and oxy images as a guide to get a better image to select the tumor ----------------------%
deoxy_guide = squeeze( msot(:,:,:,1));
deoxy_guide_max = max(deoxy_guide, [], 3);

oxy_guide = squeeze( msot(:,:,:,2));
oxy_guide_max = max(oxy_guide, [], 3);

%Comparing timepoints (motion)
%figure(106); imshowpair(I(:,:,1),I(:,:,200),"Scaling","independent");set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

%------------------- Select Tumor ------------------------%
while roi_approve == 0

    if toi_exist == 0
        %TOI
        ref_anat = I(:,:,50);
        % figure(2); imshow(ref_anat,[0, M], 'Colormap',jet); title('Draw tumor'); h=drawfreehand;
        %figure(2); imagesc(ref_anat); colormap(jet); axis image off; title('Draw tumor'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');   
        figure(103); imagesc(deoxy_guide_max); colormap("turbo"); axis image off; title('Draw tumor'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
        h=drawfreehand;
        mask_toi = createMask(h);

        str = input('Are you satisfied with the Tumor ROI? Y/N [Y]: ','s');
        if str == 'Y'
            roi_approve = 1;
        end
    end
    
end


%------------------- Remove the tumor area from the image ------------------------%
msot_wo_toi = msot;
[rows,cols, timepoints, dim4] = size(msot);

for idx=1:n_msp   
    temp1 = msot_wo_toi(:,:,:,idx);
     for i = 1:timepoints
         temp2 = temp1(:,:,i);
         temp2(mask_toi) = 0;
         temp1(:,:,i) = temp2;
     end
     msot_wo_toi(:,:,:,idx) = temp1;
end

disp('size(msot_wo_toi) = '); disp(size(msot_wo_toi));

%------------------- Procedure for AIF ------------------------%
% Map1: Maximum Intensity Projection (MIP)
% Use deoxy and oxy images as a guide to find AIF
deoxy_guide = squeeze( msot_wo_toi(:,:,:,1));
deoxy_guide_max = max(deoxy_guide, [], 3);

% figure(100); 
% subplot(2,2,1); imagesc(deoxy_guide(:,:,50)); colormap(jet); axis image off; title('deoxy 50 (use as a guide)'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(2,2,2); imagesc(deoxy_guide_max); colormap(jet); axis image off; title('deoxy max (use as a guide)'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

oxy_guide = squeeze( msot_wo_toi(:,:,:,2));
oxy_guide_max = max(oxy_guide, [], 3);

% subplot(2,2,3); imagesc(oxy_guide(:,:,50)); colormap(jet); axis image off; title('oxy 50 (use as a guide)'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(2,2,4); imagesc(oxy_guide_max); colormap(jet); axis image off; title('oxy max (use as a guide)'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

oxy_deoxy_max_mult = deoxy_guide_max .* oxy_guide_max;
% figure(102); 
% subplot(2,2,1); imagesc(oxy_deoxy_max_mult); colormap(jet); axis image off; title('Total Hb'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

% Map2: Temporal Difference Map
I_wo_toi = squeeze( msot_wo_toi(:,:,:,cl)); 
TDM = sum(abs(diff(I_wo_toi, 1, 3)),3);
% subplot(2,2,2); imagesc(TDM); colormap(jet); axis image off; title('Temporal Difference Map'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

% Combining Maps
oxy_deoxy_max_mult_TDM = oxy_deoxy_max_mult .* TDM;
% subplot(2,2,3); imagesc(oxy_deoxy_max_mult_TDM); colormap(jet); axis image off; title('Total Hb + TDM'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

%Thresholding
T=0.1*(max(oxy_deoxy_max_mult_TDM(:)));
Total_Hb_TDM_BW = imbinarize(oxy_deoxy_max_mult_TDM,T);
% subplot(2,2,4); imagesc(Total_Hb_TDM_BW); colormap(jet); axis image off; title('Total Hb + TDM - After Thresholding)'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

% Map3: Area Under the Curve
I_wo_toi_basline = I_wo_toi(:,:,1:6);
I_wo_toi_basline_mean = mean(I_wo_toi_basline,3);
I_wo_toi_injection = I_wo_toi(:,:,15:20);
I_wo_toi_injection_mean = mean(I_wo_toi_injection,3);
I_wo_toi_later = I_wo_toi(:,:,100:110);
I_wo_toi_later_mean = mean(I_wo_toi_later,3);
AUC_diff_1 = abs(imsubtract(I_wo_toi_injection_mean,I_wo_toi_basline_mean));
AUC_diff_2 = abs(imsubtract(I_wo_toi_later_mean, I_wo_toi_injection_mean));
AUC_diff_mult = AUC_diff_1 .* AUC_diff_2;
T=0.3*(max(AUC_diff_1(:)));
AUC_diff_1_BW = imbinarize(AUC_diff_1,T);
T=0.3*(max(AUC_diff_2(:)));
AUC_diff_2_BW = imbinarize(AUC_diff_2,T);
T=0.15*(max(AUC_diff_mult(:)));
AUC_diff_mult_BW = imbinarize(AUC_diff_mult,T);
% figure(103); 
% subplot(3,3,1); imagesc(I_wo_toi_basline_mean); colormap(jet); axis image off; title('I wo toi basline mean'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,2); imagesc(I_wo_toi_injection_mean); colormap(jet); axis image off; title('I wo toi injection mean'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,3); imagesc(I_wo_toi_later_mean); colormap(jet); axis image off; title('I wo toi later mean'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,4); imagesc(AUC_diff_1); colormap(jet); axis image off; title('AUC Diff 1'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,5); imagesc(AUC_diff_2); colormap(jet); axis image off; title('AUC Diff 2'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,6); imagesc(AUC_diff_mult); colormap(jet); axis image off; title('AUC Diff Mult'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,7); imagesc(AUC_diff_1_BW); colormap(jet); axis image off; title('AUC Diff 1 BW'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,8); imagesc(AUC_diff_2_BW); colormap(jet); axis image off; title('AUC Diff 2 BW'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
% subplot(3,3,9); imagesc(AUC_diff_mult_BW); colormap(jet); axis image off; title('AUC Diff Mult BW'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');

%------------------- Get Signals from candidate voxels ------------------------%
[rows, cols, timepoints ] = size(I);
time_min = squeeze( ts(:,:,:,3) ) / 60;
dd = 1;
dd_max = sum(AUC_diff_mult_BW(:));

voxel_row = zeros(1, dd_max);
voxel_col = zeros(1, dd_max);
signal = zeros(dd_max,timepoints);


for r = 1:rows
    for c = 1:cols
        if AUC_diff_mult_BW(r,c) == 1
            for t = 1:timepoints
                disp(dd + "/" + dd_max);
                voxel_row(dd) = r;
                voxel_col(dd) = c;
                signal(dd, t) = I(r, c, t);
            end
            dd = dd + 1;
         end
    end
end

save('signal.mat', 'signal')
save('voxel_row.mat', 'voxel_row')
save('voxel_col.mat', 'voxel_col')

%------------------- Apply DBSCAN (unsupervised learning) ------------------------%
minpts = 4;
kD = pdist2(normalize(signal),normalize(signal),'cityblock','Smallest',minpts);
% figure;
% plot(sort(kD(end,:)));
% title('k-distance graph')
% xlabel('Points sorted with 4th nearest distances')
% ylabel('4th nearest distances')
% grid

kD_sorted = sort(kD(end,:));
kD_sorted_smooth = smoothdata(kD_sorted,'gaussian');
[~,kD_index]  = max(abs(diff(kD_sorted_smooth, 2)));
% disp("kD_index: " + kD_index);
% disp("kD_sorted(kD_index): " + kD_sorted(kD_index));

epsilon = kD_sorted(kD_index);
%epsilon = 62.78;

dbscan_out = dbscan(normalize(signal),epsilon,minpts, 'Distance', 'cityblock');
unique_clusters = unique(dbscan_out);
numGroups = length(unique_clusters);

% disp("dbscan_out ="); disp(dbscan_out);
% disp("numGroups ="); disp(numGroups);
% disp("unique_clusters ="); disp(unique_clusters);

mean_clusters = zeros(numGroups, timepoints);

for ii = 1:numGroups
    mean_clusters(ii,:) = mean(signal(dbscan_out == unique_clusters(ii),:)); 
end

% disp("mean_clusters ="); disp(mean_clusters);

%------------------------------- Show Results----------------------------------%
figure(104);
B = arrayfun(@num2str, unique_clusters, 'UniformOutput', 0);
H = plot(time_min, mean_clusters', 'LineWidth',0.5); legend (B); title ("DBSCAN Clusters");
set(H, 'ButtonDownFcn', {@LineSelected, H})

% figure (105);
% hold on;
% plot(time_min, signal', 'Color', [0.9 0.9 0.9]);
% plot(time_min, mean_clusters');
% hold off;
% xlabel('Time (min)');
% ylabel('Signal ');

%---------------------- Isolate the AIF cluster --------------------------------%
% injection_segment = mean_clusters(:, 1:20);
% injection_segment_range = max(injection_segment, [], 2) - min(injection_segment, [], 2);
% [M, AIF_cluster_index] = max(injection_segment_range);
% disp("AIF Cluster Index =" + unique_clusters(AIF_cluster_index));

str = input('Choose an AIF cluster: ','s');
AIF_cluster_index = str2num(str) + 1;

%------------------------------- Show Results----------------------------------%
figure;
plot(time_min, mean_clusters(AIF_cluster_index,:)); title("AIF Cluster");

%------------------- Create AIF Mask ------------------------%
aif_voxels_index = find(dbscan_out ==  unique_clusters(AIF_cluster_index));
disp("aif_voxels_index ="); disp(aif_voxels_index);

aif_mask_dbscan = zeros(rows, cols);
%aif_mask_dbscan(voxel_col(aif_voxels_index),voxel_row(aif_voxels_index)) = 1;

for ind = 1:length(aif_voxels_index)
    disp(aif_voxels_index(ind) + " , " + voxel_row(aif_voxels_index(ind)) + " , " + voxel_col(aif_voxels_index(ind)));
    aif_mask_dbscan(voxel_row(aif_voxels_index(ind)),voxel_col(aif_voxels_index(ind))) = 1;
end

%disp(aif_mask_dbscan);
ref_anat = mean(I,3);
figure; imshowpair(ref_anat,aif_mask_dbscan, "method", "blend", "Scaling","independent"); title('Automatic AIF ROI'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
str = input('Are you satisfied with the Automatic AIF ROI? Y/N [Y]: ','s');
if str == 'Y'
    mask_aif = aif_mask_dbscan > 0;
else
    disp("Selecting Manual AIF ROI")
    %------------------- Select Manual AIF ------------------------%
    roi_approve = 0;
    while roi_approve == 0
        
        % AIF
        ref_anat = I(:,:,50);
        %ref_anat = mean(I,3);
        % figure(1); imshow(ref_anat,[0, M], 'Colormap',jet); title('Draw AIF'); h=drawfreehand;
        f1 = figure; imagesc(ref_anat); colormap("jet"); axis image off; title('Draw AIF'); set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');  
        h=drawfreehand;
        mask_aif = createMask(h); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Showing AIF curve from selected ROI
        I = squeeze( msot(:,:,:,3) );
        time_min = squeeze( ts(:,:,:,3) ) / 60;
        [rows,cols, timepoints ] = size(I);
        
        aif = zeros(1,timepoints);
        
            for i = 1:timepoints
                S = I(:,:,i);
                aif(i) = mean(S(mask_aif));
                %aif(i) = max(S(mask_aif)); 
                % my ROI is too big, I am removing the negative voxels around the AIF;
            end
        
        f2 = figure;
        plot(time_min, aif); title("AIF Curve from selected AIF ROI");
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        str = input('Are you satisfied with the AIF ROI? Y/N [Y]: ','s');
        if str == 'Y'
            roi_approve = 1;
        else
            close([f1 f2]);
        end
    end
end


%------------------- Select Muscle ------------------------%
%%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
% roi_approve = 0;
% while roi_approve == 0
% 
%     % Muscle
%     ref_anat = mean(I,3);
%     % figure(3); imshow(ref_anat,[0, M], 'Colormap', jet); title('Draw muscle'); h=drawfreehand;
%     figure(3); imagesc(ref_anat); colormap(jet); axis image off; title('Draw muscle');
%     set(gca, 'ydir', 'normal');
%     set(gca, 'xdir', 'reverse');
%     h=drawfreehand;
%     mask_muscle = createMask(h);
%     
%     str = input('Are you satisfied with the Muscle ROI? Y/N [Y]: ','s');
%     if str == 'Y'
%         roi_approve = 1;
%     end 
%    
% end
%%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here


% close all;

%Combine masks for saving purposes
%%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
% if toi_exist == 0
%     segm = 3*mask_muscle + 2*mask_aif + mask_toi;
% else
%     segm = segm + 3*mask_muscle + 2*mask_aif;
% end
%%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here

%%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
if toi_exist == 0
    segm = 2*mask_aif + mask_toi;
else
    segm = segm + 2*mask_aif;
end
%%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here

%% Directory for saving results
cd(saving_DIR);
filename=[datainfo.Name,'_DCE_MSOT_MSP'];
mkdir(filename);
cd(filename);

name = ["deoxy", "oxy", "ICG"];

%% DCE_MSOT fit 
for idx=1:n_msp
    
    I = squeeze( msot(:,:,:,idx) );
    time_min = squeeze( ts(:,:,:,idx) ) / 60;
    [rows,cols, timepoints ] = size(I);

    % Determining points with peak within the first 20 seconds
%     max_value = zeros(1,5);
%     for i = 1:5
%         S = I(:,:,i);
%         max_value(i) = max(S(mask_aif));
%     end
%     
%     max_value_all = max(max_value);
%     disp (max_value)
%     disp(max_value_all)
 


    
    %aif = mean_clusters(AIF_cluster_index,:);
    aif = zeros(1,timepoints);

    for i = 1:timepoints
        S = I(:,:,i);
        aif(i) = mean(S(mask_aif));
        %aif(i) = max(S(mask_aif)); 
        % my ROI is too big, I am removing the negative voxels around the AIF;
    end
    
    toi = zeros(1,timepoints);

    for i = 1:timepoints
        S = I(:,:,i);
        toi(i) = mean(S(mask_toi));
    end
    
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
%     muscle = zeros(1, timepoints);
%     for i = 1:timepoints
%         S = I(:,:,i);
%         muscle(i) = mean(S(mask_muscle));
%     end
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    
    %% Plot average signal for AIF and TOI within ROI
    figure;
    hold on;
    plot(time_min, aif, 'ro-', 'LineWidth',linewidthTK);
    plot(time_min, toi, 'bo-', 'LineWidth',linewidthTK);
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    %plot(time_min, muscle, 'go-', 'LineWidth',linewidthTK);
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    hold off;
    xlabel('Time (min)');
    ylabel('Signal ');
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    % legend ({'AIF','TOI', 'Muscle'},'Location','northeast');
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    
    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
    legend ({'AIF','TOI'},'Location','northeast');
    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here
    
     title(strcat('Signal for',{' '},name(idx)));
    savefig(strcat('Signal for',{' '},name(idx),'.fig'));
    
    
    %% Plot normalized signal and results from fitting
    figure;
    hold on;
    plot( time_min, normalize_signal2(aif,12)', 'ro-', 'LineWidth',linewidthTK); 
    plot( time_min, normalize_signal2(toi,12)', 'bo-', 'LineWidth',linewidthTK);
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    %plot( time_min, normalize_signal(muscle,12)', 'go-', 'LineWidth',linewidthTK);
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    hold off;
    xlabel('Time (min)');
    ylabel('Signal ');
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    % legend ({'AIF','TOI', 'Muscle'},'Location','northeast');
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    
    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
    legend ({'AIF','TOI'},'Location','northeast');
    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here
    title(strcat('Normalized signal for',{' '},name(idx)));
    savefig(strcat('Normalized signal for',{' '},name(idx),'.fig'));
    
    norm_signal_aif = normalize_signal2(aif,12)';
    norm_signal_toi = normalize_signal2(toi,12)';
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    %norm_signal_muscle = normalize_signal(muscle,12)';

    %save('Signal.mat', 'time_min', 'aif', 'toi', 'muscle', 'norm_signal_aif', 'norm_signal_toi', 'norm_signal_muscle');
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here

    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
    save('Signal.mat', 'time_min', 'aif', 'toi','norm_signal_aif', 'norm_signal_toi');
    %%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here

    %% Plot smoothed normalized signal and results from fitting
%     figure;
%     hold on;
%     plot( time_min, normalize_signal2(smooth_signal(aif),12)', 'ro-', 'LineWidth',linewidthTK); 
%     plot( time_min, normalize_signal2(smooth_signal(toi),12)', 'bo-', 'LineWidth',linewidthTK);
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
%     %plot( time_min, smooth_signal(muscle)', 'go-', 'LineWidth',linewidthTK);
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
%     hold off;
%     xlabel('Time (min)');
%     ylabel('Signal ');
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
%     % legend ({'AIF','TOI', 'Muscle'},'Location','northeast');
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
%     
%     %%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
%     legend ({'AIF','TOI'},'Location','northeast');
%     %%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here
%     title(strcat('Smooth Normalized signal for',{' '},name(idx)));
%     savefig(strcat('Smooth Normalized signal for',{' '},name(idx),'.fig'));
%     
%     smooth_norm_signal_aif = normalize_signal2(smooth_signal(aif),12)';
%     smooth_norm_signal_toi = normalize_signal2(smooth_signal(toi),12)';
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
%     %norm_signal_muscle = smooth_signal(muscle)';
% 
%     %save('Signal.mat', 'time_min', 'aif', 'toi', 'muscle', 'smooth_norm_signal_aif', 'smooth_norm_signal_toi', 'smooth_norm_signal_muscle');
%     %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
% 
%     %%%%%%%%%%%%%%%%%%%%% Replaced with this - Starting here
%     save('Signal.mat', 'time_min', 'aif', 'toi','smooth_norm_signal_aif', 'smooth_norm_signal_toi');
%     %%%%%%%%%%%%%%%%%%%%% Replaced with this - Ending here
%     % close all
    
    %% Fitting of average signal in ROI
    [kapp_avg_signal, kep_avg_signal, ~, pars_avg_signal] = fit_dcemsot(toi,aif,time_min);    
    name_file = strcat('Average_signal_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_avg_signal', 'kep_avg_signal','pars_avg_signal'); 

    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    % [kapp_avg_signal_muscle, kep_avg_signal_muscle,~, pars_avg_signal_muscle] = fit_dcemsot(muscle,aif,time_min);
    %name_file = strcat('Muscle_Average_signal_Wavelength_',name(idx),'.mat');
    %save(name_file, 'kapp_avg_signal_muscle', 'kep_avg_signal_muscle', 'pars_avg_signal_muscle'); 
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Ending here
    
    %% Voxel-based fit
    signal = reshape( I, [rows*cols,timepoints]);
    segm = reshape(segm, [rows*cols,1]);
    
    % Allocate space for maps
    kapp = zeros((rows*cols),1);
    kep = zeros((rows*cols),1);
    pars = zeros((rows*cols),map_no);
    
    %% Create Parameter maps using Julio function
    % TOI
    for itr = 1:rows*cols
        if segm (itr)==1
            [kapp(itr,1), kep(itr,1), ~, pars(itr,:)] = fit_dcemsot(double(signal(itr,:)),double(aif), time_min');
        end
    end
    
    kapp_nz = mean(nonzeros(kapp));
    kep_nz = mean(nonzeros(kep));
    kapp = reshape(kapp, [rows,cols]);
    kep = reshape(kep, [rows,cols]);
    pars = reshape(pars, [rows,cols,map_no]);
    segm = reshape(segm,[rows,cols]);
    rsquare_fit = pars(:,:,4);
    kapp_avg = mean(kapp(segm==1));
    kep_avg = mean(kep(segm==1));
    
    name_file = strcat('Maps_Wavelength_',name(idx),'.mat');
    save(name_file, 'kapp_nz', 'kep_nz', 'kapp_avg', 'kep_avg', 'kapp', 'kep', 'rsquare_fit', 'segm','ref_anat'); 
    
    figure;
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, kapp, 'alphadata', segm==1);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('kapp');
    colorbar;
    savefig('kapp.fig');
    % saveas(gcf, 'kapp.tiff');
    print(gcf, '-dtiff', 'kapp.tiff');
    
    figure;
    ax1 = axes;
    imagesc(ref_anat)
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax1, 'gray');
    ax2 = axes;
    imagesc(ax2, kep, 'alphadata', segm==1);
    axis image off; set(gca, 'ydir', 'normal'); set(gca, 'xdir', 'reverse');
    colormap(ax2,'jet');
    caxis(ax2, [0 1]);
    ax2.Visible = 'off';
    linkprop([ax1 ax2],'Position');
    title('kep');
    colorbar;
    savefig('kep.fig');
    % saveas(gcf, 'kep.tiff'); 
    print(gcf, '-dtiff', 'kep.tiff');
   
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    % Muscle
%     kapp = zeros((rows*cols),1);
%     kep = zeros((rows*cols),1);
%     pars = zeros((rows*cols),map_no);
%     
%     for itr = 1:rows*cols
%         if segm (itr)==3
%             [kapp(itr,1), kep(itr,1), ~, pars(itr,:)] = fit_dcemsot(double(signal(itr,:)),double(aif), time_min');
%         end
%     end
%     
%     kapp_nz = mean(nonzeros(kapp));
%     kep_nz = mean(nonzeros(kep));
%     kapp = reshape(kapp, [rows,cols]);
%     kep = reshape(kep, [rows,cols]);
%     pars = reshape(pars, [rows,cols,map_no]);
%     segm = reshape(segm,[rows,cols]);
%     rsquare_fit = pars(:,:,4);
%     kapp_avg = mean(kapp(segm==3));
%     kep_avg = mean(kep(segm==3));
%     
%     name_file = strcat('Muscle_Maps_Wavelength_',name(idx),'.mat');
%     save(name_file, 'kapp_nz', 'kep_nz', 'kapp_avg', 'kep_avg', 'kapp', 'kep', 'rsquare_fit', 'segm','ref_anat'); 
    %%%%%%%%%%%%%%%%%%%%% Commenting out muscle - by Alia Khaled on 3/14/2023 - Starting here
    
    close all;
end
end

%% Functions
function norm_sig = normalize_signal(mysignal, baseline_points)

norm_sig = mysignal - mean(mysignal(1:baseline_points));
%norm_sig = norm_sig / max(mysignal); 
norm_sig = norm_sig / max(norm_sig);
end

function norm_sig = normalize_signal2(mysignal, baseline_points)
  %mysignal_new = detrend(mysignal);
  mysignal_new = mysignal;
  norm_sig = mysignal_new - mean(mysignal_new(1:baseline_points));
  norm_sig = norm_sig / max (norm_sig);
end

function smooth_sig = smooth_signal(mysignal)
  smooth_sig = smoothdata(mysignal,'gaussian');
end

function [kapp, kep, X, pars] = fit_dcemsot(TOI,AIF,TIME)

aif_norm = normalize_signal2(AIF,12)';
toi_norm = normalize_signal2(TOI,12)';

X = zeros(length(aif_norm), 2);
X(:,1) = cumtrapz(TIME, aif_norm);
X(:,2) = -1 .* cumtrapz(TIME, toi_norm);
X(:,3) = ones(size(toi_norm));

y = toi_norm;

options = optimset('Display','Notify');
pars = lsqnonneg(X, y, options);

kapp = pars(1);
kep  = pars(2);

% Calculate predicted curve for the tumor
prediction=X*pars; 

% Calculate the Rsquare for the predicted curve and allocate to pars
pars(4)=rsquare(y,prediction);

end

function [r2, rmse] = rsquare(y,f,varargin)
% Compute coefficient of determination of data fit model and RMSE
%
% [r2 rmse] = rsquare(y,f)
% [r2 rmse] = rsquare(y,f,c)
%
% RSQUARE computes the coefficient of determination (R-square) value from
% actual data Y and model data F. The code uses a general version of 
% R-square, based on comparing the variability of the estimation errors 
% with the variability of the original values. RSQUARE also outputs the
% root mean squared error (RMSE) for the user's convenience.
%
% Note: RSQUARE ignores comparisons involving NaN values.
% 
% INPUTS
%   Y       : Actual data
%   F       : Model fit
%
% OPTION
%   C       : Constant term in model
%             R-square may be a questionable measure of fit when no
%             constant term is included in the model.
%   [DEFAULT] TRUE : Use traditional R-square computation
%            FALSE : Uses alternate R-square computation for model
%                    without constant term [R2 = 1 - NORM(Y-F)/NORM(Y)]
%
% OUTPUT 
%   R2      : Coefficient of determination
%   RMSE    : Root mean squared error
%
% EXAMPLE
%   x = 0:0.1:10;
%   y = 2.*x + 1 + randn(size(x));
%   p = polyfit(x,y,1);
%   f = polyval(p,x);
%   [r2 rmse] = rsquare(y,f);
%   figure; plot(x,y,'b-');
%   hold on; plot(x,f,'r-');
%   title(strcat(['R2 = ' num2str(r2) '; RMSE = ' num2str(rmse)]))
%   
% Jered R Wells
% 11/17/11
% jered [dot] wells [at] duke [dot] edu
%
% v1.2 (02/14/2012)
%
% Thanks to John D'Errico for useful comments and insight which has helped
% to improve this code. His code POLYFITN was consulted in the inclusion of
% the C-option (REF. File ID: #34765).

if isempty(varargin); c = true; 
elseif length(varargin)>1; error 'Too many input arguments';
elseif ~islogical(varargin{1}); error 'C must be logical (TRUE||FALSE)'
else c = varargin{1}; 
end

% Compare inputs
if ~all(size(y)==size(f)); error 'Y and F must be the same size'; end

% Check for NaN
tmp = ~or(isnan(y),isnan(f));
y = y(tmp);
f = f(tmp);

if c; r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
else r2 = 1 - sum((y(:)-f(:)).^2)/sum((y(:)).^2);
    if r2<0
    % http://web.maths.unsw.edu.au/~adelle/Garvan/Assays/GoodnessOfFit.html
        warning('Consider adding a constant term to your model') %#ok<WNTAG>
        r2 = 0;
    end
end

rmse = sqrt(mean((y(:) - f(:)).^2));

end

function LineSelected(ObjectH, EventData, H)
set(ObjectH, 'LineWidth', 2.5);
set(H(H ~= ObjectH), 'LineWidth', 0.5);
end
